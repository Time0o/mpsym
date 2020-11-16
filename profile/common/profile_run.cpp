#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "util.hpp"

#include "profile_run.hpp"
#include "profile_timer.hpp"

namespace
{

class TmpFile
{
public:
  TmpFile(std::string const &content, std::string const &fname = "")
  : _fname(fname)
  {
    if (_fname.empty()) {
      char tmpname[6] = {'X', 'X', 'X', 'X', 'X', 'X'};
      if (mkstemp(tmpname) == -1)
        throw std::runtime_error("failed to create temporary file");

      _fname = tmpname;
    }

    _fname += ".g";

    _f = std::ofstream(name());

    if (_f.fail())
      throw std::runtime_error("failed to create temporary file");

    _f << content;
    _f.flush();
  }

  ~TmpFile()
  {
    if (_valid)
      std::remove(name());
  }

  TmpFile(TmpFile const &) = delete;
  TmpFile & operator=(TmpFile const &) = delete;

  TmpFile(TmpFile &&other)
  : _f(std::move(other._f)),
    _fname(std::move(other._fname))
  { other._valid = false; }

  TmpFile & operator=(TmpFile &&) = delete;

  char const *name()
  { return _fname.c_str(); }

private:
  std::ofstream _f;
  std::string _fname;
  bool _valid = true;
};

std::string common_functions()
{
  return R"(
GraphAutoms:=function(edges, partition, n)
  return AutGroupGraph(EdgeOrbitsGraph(Group(()), edges, n), partition);
end;

ReduceGroup:=function(G, n)
  local gens_, gens;

  if IsTrivial(G) then
    return G;
  fi;

  gens_:=GeneratorsOfGroup(G);
  gens:=ShallowCopy(gens_);
  Apply(gens, function(g) return RestrictedPerm(g, [1..n]); end);
  return Group(gens);
end;

FixedPointWreathProduct:=function(G, nG, H, nH)
  if LargestMovedPoint(G) <> nG or LargestMovedPoint(H) <> nH then
    Error("TODO: consider fixed points");
  fi;

  return WreathProduct(G, H);
end;
)" + 1;
}

std::string build_load_script(
  std::initializer_list<std::string> packages,
  std::initializer_list<profile::Preload> preloads)
{
  std::stringstream ss;

  for (auto const &package : packages)
    ss << "LoadPackage(\"" << package << "\");\n";

  ss << '\n' << common_functions() << '\n';

  for (auto const &preload : preloads) {
    std::string preload_file(std::get<0>(preload));
    bool preload_compile = std::get<2>(preload);

    if (preload_compile)
      ss << "LoadDynamicModule(\"./" << preload_file << ".la.so" << "\");\n";
    else
      ss << "Read(\"" << preload_file << ".g" << "\");\n";
  }

  return ss.str();
}

std::string build_script(
  std::initializer_list<std::string> packages,
  std::initializer_list<profile::Preload> preloads,
  std::string const &script,
  unsigned num_discarded_runs,
  unsigned num_runs)
{
  std::stringstream ss;

  ss << build_load_script(packages, preloads) << '\n';

  ss << "_ts:=[];\n";
  ss << "for _r in [1.." << num_discarded_runs + num_runs << "] do\n";
  ss << "  _start:=NanosecondsSinceEpoch();\n";
  ss << script;
  ss << "  if _r > " << num_discarded_runs << " then\n";
  ss << "    Add(_ts, NanosecondsSinceEpoch() - _start);\n";
  ss << "  fi;\n";
  ss << "od;\n";
  ss << "Print(\"RESULT: \", _ts, \"\\n\");\n";
  ss << "Print(\"END\\n\");\n";

  return ss.str();
}

std::string build_wrapper_script(
  std::initializer_list<std::string> packages,
  std::initializer_list<profile::Preload> preloads,
  std::string lib)
{
  std::stringstream ss;

  ss << build_load_script(packages, preloads) << '\n';

  ss << "LoadDynamicModule(\"./" << lib << ".la.so" << "\");";

  return ss.str();
}

void dup_fd(int from, int to)
{
  for (;;) {
    if (dup2(from, to) == -1) {
      if (errno == EINTR)
        continue;
      else
        throw std::runtime_error("dup failed");
    } else {
      break;
    }
  }
}

void connect_stream(int stream, int *pipe)
{
  dup_fd(pipe[1], stream);
  close(pipe[1]);
  close(pipe[0]);
}

void redirect_stream(int stream, int to)
{ dup_fd(to, stream); }

std::string read_output(int from, bool echo)
{
  static char buf[256];

  std::string res;

  for (;;) {
    auto count = read(from, buf, sizeof(buf));

    if (count == -1) {
      if (errno == EINTR)
        continue;
      else
        throw std::runtime_error("read failed");
      break;
    } else if (count == 0) {
      break;
    } else {
      std::string block(buf, buf + count);

      if (block.find("END") != std::string::npos) {
        res += block.substr(0u, block.size() - std::strlen("END") - 1u);
        break;
      }

      res += block;

      if (block.find("RESULT") != std::string::npos)
        echo = false;

      if (echo) {
        auto echo_block(block);

        echo_block.erase(std::remove(echo_block.begin(), echo_block.end(), ';'),
                         echo_block.end());

        std::cout << echo_block << std::flush;
      }
    }
  }

  return res;
}

template<typename FUNC>
std::string run_in_child(FUNC &&f,
                         int *output_pipe = nullptr,
                         bool hide_stdout = false,
                         bool hide_stderr = false)
{
  // fork child process
  pid_t child;
  switch ((child = fork())) {
  case -1:
    throw std::runtime_error("failed to fork child process");
  case 0:
    {
      // connect stdout/stderr pipes
      if (output_pipe)
        connect_stream(STDOUT_FILENO, output_pipe);

      // hide stdout/stderr
      int dev_null = open("/dev/null", O_WRONLY);

      if (dev_null != -1) {
        if (!output_pipe && hide_stdout)
          redirect_stream(STDOUT_FILENO, dev_null);

        if (hide_stderr)
          redirect_stream(STDERR_FILENO, dev_null);
      }

      close(dev_null);

      // run function in child process
      if (f() == -1)
        _Exit(EXIT_FAILURE);
    }
  }

  // read output
  std::string output;

  if (output_pipe) {
    output = read_output(output_pipe[0], !hide_stdout);

    if (!hide_stdout)
      std::cout << std::endl;

    close(output_pipe[1]);
    close(output_pipe[0]);
  } else {
    output = "";
  }

  // check child process exit status
  int status;
  if (waitpid(child, &status, 0) == -1)
    throw std::runtime_error("waitpid failed");

  if (!WIFEXITED(status))
    throw std::runtime_error("child process did not terminate normally");

  // return output
  return output;
}

std::string compile_script(std::string script, std::string file)
{
  TmpFile f_tmp(script, file);

  run_in_child([&]{
      return execlp("gac", "gac", "-d", f_tmp.name(), nullptr);
    },
    nullptr,
    true,
    true);

  return file + ".la.so";
}

void remove_libs(std::vector<std::string> const &libs)
{
  for (auto const &lib : libs)
    std::remove(lib.c_str());
}

std::string clean_output(std::string const &output)
{
  std::size_t i = 0u, j = output.size() - 1u;

  while (std::isspace(output[i]))
    ++i;

  while (std::isspace(output[j]))
    --j;

  return output.substr(i, j - i + 1u);
}

std::string compress_output(std::string const &output)
{
  auto res(output);

  for (char space : " \n\\")
    res.erase(std::remove(res.begin(), res.end(), space), res.end());

  return res;
}

std::vector<std::string> parse_output(std::string const &output_,
                                      unsigned num_runs,
                                      std::vector<double> *ts)
{
  using mpsym::util::split;
  using mpsym::util::stof;

  auto output(compress_output(clean_output(output_)));

  // parse output
  auto output_split(split(output, ";"));

  std::vector<std::string> output_vals(
    output_split.begin(),
    output_split.begin() + (output_split.size() - 1u) / num_runs);

  if (ts) {
    auto parse_output_times = [](std::string output_times){
      output_times = output_times.substr(std::strlen("RESULT:"));
      output_times = output_times.substr(1u, output_times.size() - 2u);

      return split(output_times, ",");
    };

    for (auto const &output_time : parse_output_times(output_split.back()))
      ts->push_back(stof<double>(output_time) / 1e9);
  }

  return output_vals;
}

}

namespace profile
{

std::vector<std::string> run_gap(
  std::initializer_list<std::string> packages,
  std::initializer_list<Preload> preloads,
  std::string const &script,
  unsigned num_discarded_runs,
  unsigned num_runs,
  bool hide_output,
  bool hide_errors,
  bool compile,
  std::vector<double> *ts)
{
  // list of shared libs to be removed later
  std::vector<std::string> libs;

  // create preload files
  std::vector<TmpFile> f_preloads;

  for (auto const &preload : preloads) {
    std::string preload_file(std::get<0>(preload));
    std::string preload_content(std::get<1>(preload));
    bool preload_compile = std::get<2>(preload);

    if (preload_compile)
      libs.push_back(compile_script(preload_content, preload_file));
    else
      f_preloads.emplace_back(preload_content, preload_file);
  }

  // create gap script
  std::string script_main;

  if (compile) {
    // compile script
    auto script_compiled(build_script({},
                                      {},
                                      script,
                                      num_discarded_runs,
                                      num_runs));

    libs.push_back(compile_script(script_compiled, "compiled"));

    // construct main script loading compiled shared object
    script_main = build_wrapper_script(packages, preloads, "compiled");

  } else {
    script_main = build_script(packages,
                               preloads,
                               script,
                               num_discarded_runs,
                               num_runs);
  }

  TmpFile f_script(script_main, "script");

  // create pipe for capturing gaps output
  int output_pipe[2];
  if (pipe(output_pipe) == -1)
    throw std::runtime_error("failed to create pipe");

  // run gap in a child process
  auto output(run_in_child([&]{
      return execlp("gap", "gap", "--nointeract", "-q", f_script.name(), nullptr);
    },
    output_pipe,
    hide_output,
    hide_errors));

  // remove compiled shared objects
  if (compile)
    remove_libs(libs);

  // parse and return output
  return parse_output(output, num_runs, ts);
}

void dump_runs(std::vector<double> const &ts, bool summarize)
{
  if (summarize) {
    double t_mean, t_stddev;
    mpsym::util::mean_stddev(ts, &t_mean, &t_stddev);

    result("Mean:", t_mean, "s");
    result("Stddev:", t_stddev, "s");

  } else {
    result("Runtimes:");
    for (double t : ts)
      result(t, "s");
  }
}

} // namespace profile
