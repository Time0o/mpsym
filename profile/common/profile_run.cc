#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include "profile_run.h"
#include "profile_timer.h"

namespace
{

std::string build_script(std::initializer_list<std::string> packages,
                         std::string const &script,
                         unsigned num_discarded_runs,
                         unsigned num_runs)
{
  std::stringstream ss;

  for (auto const &package : packages)
    ss << "LoadPackage(\"" << package << "\");\n";

  ss << "_ts:=[];\n";
  ss << "for _r in [1.." << num_discarded_runs + num_runs << "] do\n";
  ss << "  _start:=NanosecondsSinceEpoch();\n";
  ss << script;
  ss << "  if _r > " << num_discarded_runs << " then\n";
  ss << "    Add(_ts, NanosecondsSinceEpoch() - _start);\n";
  ss << "  fi;\n";
  ss << "od;\n";
  ss << "Print(_ts, \"\\n\");\n";

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

      if (echo)
        std::cout << block << std::flush;

      res += block;
    }
  }

  return res;
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

}

namespace profile
{

std::vector<std::string> run_gap(std::initializer_list<std::string> packages,
                                 std::string const &script,
                                 unsigned num_discarded_runs,
                                 unsigned num_runs,
                                 bool hide_output,
                                 bool hide_errors,
                                 std::vector<double> *ts)
{
  // create temporary gap script

  char ftmp[] = {'X', 'X', 'X', 'X', 'X', 'X'};
  if (mkstemp(ftmp) == -1)
    throw std::runtime_error("failed to create temporary file");

  class FileRemover
  {
  public:
    FileRemover(std::string const &f) : _f(f) {}
    ~FileRemover() { std::remove(_f.c_str()); }
  private:
    std::string _f;
  } file_remover(ftmp);

  std::ofstream f(ftmp);
  if (f.fail())
    throw std::runtime_error("failed to create temporary file");

  f << build_script(packages, script, num_discarded_runs, num_runs);

  f.flush();

  // create pipe for capturing gaps output

  int fds[2];
  if (pipe(fds) == -1)
    throw std::runtime_error("failed to create pipe");

  // run gap in a child process

  timer_start();

  pid_t child;
  switch ((child = fork())) {
  case -1:
    throw std::runtime_error("failed to fork child process");
  case 0:
    {
      dup_fd(fds[1], STDOUT_FILENO);

      fcntl(STDOUT_FILENO, F_SETFL, fcntl(STDOUT_FILENO, F_GETFL) | O_DIRECT);

      close(fds[1]);
      close(fds[0]);

      if (hide_errors) {
        int dev_null = open("/dev/null", O_WRONLY);
        if (dev_null != -1)
          dup_fd(dev_null, STDERR_FILENO);
      }

      if (execlp("gap", "gap", "--nointeract", "-q", ftmp, nullptr) == -1)
        _Exit(EXIT_FAILURE);

      _Exit(EXIT_SUCCESS);
    }
  }

  close(fds[1]);

  // get output

  auto output(compress_output(clean_output(read_output(fds[0], !hide_output))));

  // parse output
  auto output_split(split(output, ";"));

  std::vector<std::string> output_vals(
    output_split.begin(),
    output_split.begin() + (output_split.size() - 1u) / num_runs);

  if (ts) {
    auto parse_output_times = [](std::string const &output_times){
      return split(output_times.substr(1u, output_times.size() - 2u), ",");
    };

    for (auto const &output_time : parse_output_times(output_split.back()))
      ts->push_back(stox<double>(output_time) / 1e9);
  }

  close(fds[0]);

  return output_vals;
}

} // namespace profile
