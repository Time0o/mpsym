#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <vector>

#include <getopt.h>
#include <libgen.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "perm.h"
#include "perm_group.h"
#include "util.h"

namespace boost {
  // permlib needs this to compile and I couldn't be bothered to build this
  // program against an appropriate version of boost instead
  template<typename ForwardIt>
  ForwardIt next(
    ForwardIt it,
    typename std::iterator_traits<ForwardIt>::difference_type n = 1) {

    return std::next(it, n);
  }
}

#include <permlib/permlib_api.h>

#include <permlib/construct/random_schreier_sims_construction.h>
#include <permlib/generator/bsgs_random_generator.h>
#include <permlib/transversal/explicit_transversal.h>
#include <permlib/transversal/shallow_schreier_tree_transversal.h>


namespace {

std::string progname;

// output

void usage(std::ostream &s)
{
  char const *opts[] = {
    "[-h|--help]",
    "-i|--implementation        {mpsym|permlib|gap}",
    "-s|--schreier-sims         {deterministic|random}",
    "[-t|--transversal-storage] {explicit|schreier-tree|random-schreier-tree}",
    "[-c|--num-cycles]",
    "[-r|--num-runs]",
    "[-v|--verbose]",
    "GROUPS"
  };

  s << "usage: " << progname << '\n';
  for (char const *opt : opts)
    s << "  " << opt << '\n';
}

template<typename Arg, typename... Args>
void error(Arg &&arg, Args&&... args) {
  std::cerr << progname << ": error: ";

  std::cerr << std::forward<Arg>(arg);

  using expander = int[];
  (void)expander{0, (void(std::cerr << ' ' << std::forward<Args>(args)), 0)...};

  std::cerr << '\n';
}

// argument parsing

enum LibraryImpl {
  LIBRARY_UNSET,
  LIBRARY_MPSYM,
  LIBRARY_PERMLIB,
  LIBRARY_GAP
};

enum SchreierSimsImpl {
  SCHREIER_SIMS_UNSET,
  SCHREIER_SIMS_DETERMINISTIC,
  SCHREIER_SIMS_RANDOM
};

enum TransversalImpl {
  TRANSVERSAL_UNSET,
  TRANSVERSAL_EXPLICIT,
  TRANSVERSAL_SCHREIER_TREE,
  TRANSVERSAL_SHALLOW_SCHREIER_TREE
};

std::unordered_map<std::string, LibraryImpl> library_impls {
  { "mpsym", LIBRARY_MPSYM },
  { "permlib", LIBRARY_PERMLIB },
  { "gap", LIBRARY_GAP }
};

std::unordered_map<std::string, SchreierSimsImpl> schreier_sims_impls {
  { "deterministic", SCHREIER_SIMS_DETERMINISTIC },
  { "random", SCHREIER_SIMS_RANDOM }
};

std::unordered_map<std::string, TransversalImpl> transversal_impls {
  { "explicit", TRANSVERSAL_EXPLICIT },
  { "schreier-tree", TRANSVERSAL_SCHREIER_TREE },
  { "shallow-schreier-tree", TRANSVERSAL_SHALLOW_SCHREIER_TREE }
};

template<typename T>
T choose_impl(std::unordered_map<std::string, T> params,
              std::string const &choice)
{
  auto it = params.find(choice);
  if (it == params.end())
    throw std::invalid_argument("invalid parameter choice");

  return it->second;
}

void stoi_strict(std::string const &str, int *i)
{
  bool success = true;

  std::size_t idx;
  try {
    *i = std::stoi(str, &idx);
    success = idx == str.size();
  } catch (...) {
    success = false;
  }

  if (!success)
    throw std::invalid_argument("stoi failed");
}


// group construction

void make_mpsym_perm_group(SchreierSimsImpl schreier_sims_impl,
                           TransversalImpl transversal_impl,
                           unsigned degree,
                           std::vector<cgtl::Perm> const &gens,
                           int num_cycles)
{
  using cgtl::PermGroup;

  PermGroup::ConstructionMethod constr =
    schreier_sims_impl == SCHREIER_SIMS_DETERMINISTIC ?
     PermGroup::SCHREIER_SIMS :
    schreier_sims_impl == SCHREIER_SIMS_RANDOM ?
      PermGroup::SCHREIER_SIMS_RANDOM :
    throw std::logic_error("unreachable");

  PermGroup::TransversalStorageMethod transv =
    transversal_impl == TRANSVERSAL_EXPLICIT ?
      PermGroup::EXPLICIT_TRANSVERSALS :
    transversal_impl == TRANSVERSAL_SCHREIER_TREE ?
      PermGroup::SCHREIER_TREES :
    transversal_impl == TRANSVERSAL_SHALLOW_SCHREIER_TREE ?
      PermGroup::SHALLOW_SCHREIER_TREES :
    throw std::logic_error("unreachable");

  for (int i = 0; i < num_cycles; ++i)
    PermGroup g(degree, gens, constr, transv);
}

template <typename T>
struct TypeTag { using type = T; };

std::variant<
  TypeTag<permlib::ExplicitTransversal<permlib::Permutation>>,
  TypeTag<permlib::SchreierTreeTransversal<permlib::Permutation>>,
  TypeTag<permlib::ShallowSchreierTreeTransversal<permlib::Permutation>>
>
permlib_transversal_type(TransversalImpl transversal_impl)
{
  using namespace permlib;

  switch (transversal_impl) {
    case TRANSVERSAL_EXPLICIT:
      return TypeTag<ExplicitTransversal<Permutation>>{};
    case TRANSVERSAL_SCHREIER_TREE:
      return TypeTag<SchreierTreeTransversal<Permutation>>{};
    case TRANSVERSAL_SHALLOW_SCHREIER_TREE:
      return TypeTag<ShallowSchreierTreeTransversal<Permutation>>{};
    default:
      throw std::logic_error("unreachable");
  }
}

void make_permlib_perm_group(
  SchreierSimsImpl schreier_sims_impl,
  TransversalImpl transversal_impl,
  unsigned degree,
  std::vector<permlib::Permutation::ptr> const &gens,
  int num_cycles)
{
  using namespace permlib;

  std::visit([&](auto transv_type_) {
    using transv_type = typename decltype(transv_type_)::type;

    if (schreier_sims_impl == SCHREIER_SIMS_DETERMINISTIC) {
      SchreierSimsConstruction<Permutation, transv_type>
      construction(degree);

      for (int i = 0; i < num_cycles; ++i)
        construction.construct(gens.begin(), gens.end());

    } else { // SCHREIER_SIMS_RANDOM
      BSGS<Permutation, transv_type> bsgs(degree);

      std::unique_ptr<BSGSRandomGenerator<Permutation, transv_type>>
      random_generator(new BSGSRandomGenerator<Permutation, transv_type>(bsgs));

      RandomSchreierSimsConstruction<Permutation, transv_type>
      construction(degree, random_generator.get());

      bool guaranteed = true;
      for (int i = 0; i < num_cycles; ++i)
        construction.construct(gens.begin(), gens.end(), guaranteed);
    }
  }, permlib_transversal_type(transversal_impl));
}

// cpu timer

double time_child_acc;

void time_child_start() {}

bool time_child_stop(pid_t child, double *t)
{
  int status;
  waitpid(child, &status, 0);

  if (!WIFEXITED(status) || WEXITSTATUS(status) != EXIT_SUCCESS)
    return false;

  struct tms tms;
  times(&tms);

  *t = static_cast<double>(tms.tms_cutime) / sysconf(_SC_CLK_TCK) - time_child_acc;
  time_child_acc += *t;

  return true;
}


// realtime timer

std::chrono::high_resolution_clock::time_point time_realtime_begin;

void time_realtime_start()
{
  time_realtime_begin = std::chrono::high_resolution_clock::now();
}

void time_realtime_stop(double *t)
{
  auto delta = std::chrono::high_resolution_clock::now() - time_realtime_begin;
  auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(delta);

  *t = static_cast<double>(ns.count()) / 10e9;
}

// timer wrapper functions

bool time_realtime_enabled = false;

pid_t timer_start() {
  if (time_realtime_enabled) {
    time_realtime_start();
    return 0;
  } else {
    time_child_start();
    return fork();
  }
}

double timer_stop(pid_t maybe_child) {
  double t;
  if (time_realtime_enabled) {
    time_realtime_stop(&t);
  } else {
    if (!time_child_stop(maybe_child, &t)) {
      error("the forked child process terminated prematurely");
      return false;
    }
  }

  return t;
}

template<LibraryImpl L>
bool run_cpp(SchreierSimsImpl schreier_sims_impl,
             TransversalImpl transversal_impl,
             std::string const &generators,
             int num_cycles,
             double *t)
{
  // validate generator expression
  std::string perm = R"((\(\)|(\((\d+,)+\d+\))+))";

  std::regex re("\\[(" + perm + ",)*(" + perm + ")?\\]");

  if (!std::regex_match(generators, re)) {
    error("malformed generator expression");
    return false;
  }

  // extract permutation expressions
  std::vector<std::string> gen_strs;

  std::size_t pos = 1, nextpos;
  for (;;) {
    nextpos = generators.find("),", pos);

    if (nextpos == std::string::npos) {
      gen_strs.push_back(generators.substr(pos, generators.size() - pos - 1));
      break;
    } else {
      gen_strs.push_back(generators.substr(pos, nextpos - pos + 1));
    }

    pos = nextpos + 2;
  }

  // parse permutation expressions
  unsigned degree = 0;
  std::vector<std::vector<std::vector<unsigned>>> gens_;

  for (auto const &gen_str : gen_strs) {
    std::vector<unsigned> cycle;
    std::vector<std::vector<unsigned>> perm;

    int n_beg = -1;
    for (int i = 0; i < static_cast<int>(gen_str.size()); ++i) {
      char c = gen_str[i];

      switch (c) {
      case '(':
        cycle.clear();
        break;
      case ',':
      case ')':
        {
          int n = 0;

          try {
            stoi_strict(gen_str.substr(n_beg, i - n_beg), &n);
          } catch (std::invalid_argument const &) {
            assert(false);
          }

          degree = std::max(degree, static_cast<unsigned>(n));

          cycle.push_back(n);
          if (c == ')')
            perm.push_back(cycle);

          n_beg = -1;
        }
        break;
      default:
        if (n_beg == -1)
          n_beg = i;
      }
    }

    gens_.push_back(perm);
  }

  typedef typename std::conditional<
    L == LIBRARY_MPSYM, cgtl::Perm,
    permlib::Permutation::ptr
  >::type gen_type;

  std::vector<gen_type> gens(gens_.size());

  // construct generators
  if constexpr (L == LIBRARY_MPSYM) {
    for (auto i = 0u; i < gens_.size(); ++i)
      gens[i] = cgtl::Perm(degree, gens_[i]);

  } else { // LIBRARY_PERMLIB
    for (auto i = 0u; i < gens_.size(); ++i) {
      auto &gen(gens_[i]);

      std::stringstream gen_str;
      for (auto j = 0u; j < gen.size(); ++j) {
        auto &cycle(gens_[i][j]);
        for (auto k = 0u; k < gens_[i][j].size(); ++k) {
          gen_str << i << cycle[k];

          if (k == cycle.size() - 1) {
            if (j < gen.size() - 1)
              gen_str << ", ";
          } else {
            gen_str << " ";
          }
        }
      }

      gens[i] = permlib::Permutation::ptr(
        new permlib::Permutation(degree, gen_str.str()));
    }
  }

  // run group creation in child process
  pid_t maybe_child;
  switch ((maybe_child = timer_start())) {
  case -1:
    error("failed to fork child process");
    return false;
    break;
  case 0:
    {
      if constexpr (L == LIBRARY_MPSYM) {
        make_mpsym_perm_group(schreier_sims_impl,
                              transversal_impl,
                              degree,
                              gens,
                              num_cycles);
      } else { // LIBRARY_PERMLIB
        make_permlib_perm_group(schreier_sims_impl,
                                transversal_impl,
                                degree,
                                gens,
                                num_cycles);
      }

      if (!time_realtime_enabled)
        _Exit(EXIT_SUCCESS);
    }
    break;
  }

  *t = timer_stop(maybe_child);

  return true;
}

bool run_gap(std::string const &generators, int num_cycles, double *t)
{
  char ftmp[] = {'X', 'X', 'X', 'X', 'X', 'X'};
  if (mkstemp(ftmp) == -1) {
    error("failed to create temporary file");
    return false;
  }

  class FileRemover {
  public:
    FileRemover(std::string const &f) : _f(f) {}
    ~FileRemover() { std::remove(_f.c_str()); }
  private:
    std::string _f;
  } file_remover(ftmp);

  std::ofstream f(ftmp);
  if (f.fail()) {
    error("failed to create temporary file");
    return false;
  }

  f << "for i in [1.." << num_cycles << "] do\n";
  f << "  StabChain(Group(" + generators + "));\n";
  f << "od;\n";
  f.flush();

  pid_t maybe_child;
  switch ((maybe_child = timer_start())) {
  case -1:
    error("failed to fork child process");
    return false;
    break;
  case 0:
    {
      if (execlp("gap", "gap", "--nointeract", "-q", ftmp, nullptr) == -1) {
        error("failed to exec gap");
        _Exit(EXIT_FAILURE);
      }

      if (!time_realtime_enabled)
        _Exit(EXIT_SUCCESS);
    }
    break;
  }

  *t = timer_stop(maybe_child);

  return true;
}

}

int main(int argc, char **argv)
{
  progname = basename(argv[0]);

  struct option long_options[] = {
    {"help",                no_argument,       0,       'h'},
    {"implementation",      required_argument, 0,       'i'},
    {"schreier-sims",       required_argument, 0,       's'},
    {"transversal-storage", required_argument, 0,       't'},
    {"num-cyles",           required_argument, 0,       'c'},
    {"num-runs",            required_argument, 0,       'r'},
    {"realtime-clock",      no_argument,       0,        1 },
    {"verbose",             no_argument,       0,       'v'},
    {nullptr,               0,                 nullptr,  0 }
  };

  LibraryImpl library_impl = LIBRARY_UNSET;
  SchreierSimsImpl schreier_sims_impl = SCHREIER_SIMS_UNSET;
  TransversalImpl transversal_impl = TRANSVERSAL_UNSET;

  std::string groups;

  int num_cycles = 1;
  int num_runs = 1;
  bool verbose = false;

  // TODO: detect superfluous non-argument options

  for (;;) {
    int c = getopt_long(argc, argv, "hi:s:t:r:c:v", long_options, nullptr);
    if (c == -1)
      break;

    try {
      switch(c) {
      case 'h':
        usage(std::cout);
        return EXIT_SUCCESS;
      case 'i':
        library_impl = choose_impl(library_impls, optarg);
        break;
      case 's':
        schreier_sims_impl = choose_impl(schreier_sims_impls, optarg);
        break;
      case 't':
        transversal_impl = choose_impl(transversal_impls, optarg);
        break;
      case 'c':
        stoi_strict(optarg, &num_cycles);
        break;
      case 'r':
        stoi_strict(optarg, &num_runs);
        break;
      case 'v':
        verbose = true;
        break;
      case 1:
        time_realtime_enabled = true;
        break;
      default:
        return EXIT_FAILURE;
      }
    } catch (std::invalid_argument const &) {
      error("invalid argument to", long_options[optind].name);
      return EXIT_FAILURE;
    }
  }

  if (optind == argc) {
    usage(std::cerr);
    error("GROUPS argument is mandatory");
    return EXIT_FAILURE;
  }

  groups = argv[optind];

  if (library_impl == LIBRARY_UNSET) {
    usage(std::cerr);
    error("--implementation option is mandatory");
    return EXIT_FAILURE;
  }

  if (schreier_sims_impl == SCHREIER_SIMS_UNSET) {
    usage(std::cerr);
    error("--schreier-sims option is mandatory");
    return EXIT_FAILURE;
  }

  if (library_impl != LIBRARY_GAP && transversal_impl == TRANSVERSAL_UNSET) {
    usage(std::cerr);
    error("--transversal-storage option is mandatory when not using gap");
    return EXIT_FAILURE;
  }

  std::ifstream f(groups);
  if (f.fail()) {
    error("failed to open " + groups);
    return EXIT_FAILURE;
  }

  if (verbose)
    Timer::enabled = true;

  std::regex re("degree:(\\d+),order:(\\d+),gens:(.*)");
  std::smatch m;

  std::string line;
  int lineno = 1;
  while (std::getline(f, line)) {
    if (!std::regex_match(line, m, re)) {
      error("failed to parse line no.", lineno, "in", groups);
      return EXIT_FAILURE;
    }

    auto degree = m[1];
    auto order = m[2];
    auto generators = m[3];

    if (verbose) {
      std::cout << "profiling group " << lineno
                << " with degree " << degree
                << " and order " << order
                << std::endl;
    }

    std::vector<double> ts;
    for (int r = 0; r < num_runs; ++r) {
      if (verbose)
        std::cout << "run " << r + 1 << '/' << num_runs << std::endl;

      double t;
      if (library_impl == LIBRARY_GAP) {
        if (!run_gap(generators, num_cycles, &t)) {
          error("failed to run gap");
          return EXIT_FAILURE;
        }
      } else {
        bool success = true;

        if (library_impl == LIBRARY_MPSYM) {
          success = run_cpp<LIBRARY_MPSYM>(schreier_sims_impl,
                                           transversal_impl,
                                           generators,
                                           num_cycles,
                                           &t);
        } else {
          success = run_cpp<LIBRARY_PERMLIB>(schreier_sims_impl,
                                             transversal_impl,
                                             generators,
                                             num_cycles,
                                             &t);
        }

        if (!success) {
          error("profiling failed");
          return EXIT_FAILURE;
        }
      }

      ts.push_back(t);
    }

    double t_mean, t_stddev;
    util::mean_stddev(ts, &t_mean, &t_stddev);

    std::cout.precision(3);
    std::cout << "mean: " << std::scientific << t_mean << "s, "
              << "stddev: " << std::scientific << t_stddev << "s"
              << std::endl;

    ++lineno;
  }

  if (f.bad()) {
    error("failed to read " + groups);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
