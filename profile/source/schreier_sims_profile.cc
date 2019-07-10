#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>

#include <getopt.h>
#include <libgen.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "perm.h"
#include "perm_group.h"

using cgtl::Perm;
using cgtl::PermGroup;


namespace {

std::string progname;

// output

void usage(std::ostream &s)
{
  std::string header = "usage: " + progname;
  std::string pad = std::string(header.size() + 1, ' ');

  char const *opts[] = {
    "[-h|--help]",
    "-s|--schreier-sims",
    "[-t|--transversal-storage]",
    "[-n|--num-runs]",
    "[-v|--verbose]",
    "GROUPS"
  };

  s << header << ' ' << opts[0] << '\n';
  for (auto i = 1u; i < sizeof(opts) / sizeof(opts[0]); ++i)
    s << pad << opts[i] << '\n';
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

std::vector<std::string> schreier_sims_choices {
  "standard", "random", "gap"
};

std::vector<std::string> transversal_storage_choices {
  "explicit", "schreier-trees", "shallow-schreier-trees"
};

bool valid_choice(std::vector<std::string> const &choices,
                  std::string const &opt)
{
  return std::find(choices.begin(), choices.end(), opt) != choices.end();
}

// group construction

bool time_child(pid_t child, double *t, double *t_acc)
{
  int status;
  waitpid(child, &status, 0);

  if (!WIFEXITED(status) || WEXITSTATUS(status) != EXIT_SUCCESS)
    return false;

  struct tms tms;
  times(&tms);

  *t = static_cast<double>(tms.tms_cutime) / sysconf(_SC_CLK_TCK) - *t_acc;
  *t_acc += *t;

  return true;
}

bool run_cpp(std::string const &generators,
             std::string const &schreier_sims,
             std::string const &transversal_storage,
             double *t)
{
  static double t_acc = 0.0;

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

  for (auto const &gen_str : gen_strs) {
    for (char c : gen_str) {
      if (c >= '1' && c <= '9')
        degree = std::max(degree, static_cast<unsigned>(c - '0'));
    }
  }

  std::vector<Perm> gens;

  for (auto const &gen_str : gen_strs) {
    std::vector<unsigned> cycle;
    std::vector<std::vector<unsigned>> perm;

    for (char c : gen_str) {
      switch (c) {
      case '(':
        cycle.clear();
        break;
      case ')':
        perm.push_back(cycle);
        break;
      case ',':
        continue;
      default:
        cycle.push_back(static_cast<unsigned>(c - '0'));
      }
    }

    gens.push_back(Perm(degree, perm));
  }

  pid_t child;
  switch((child = fork())) {
  case -1:
    error("failed to fork child process");
    return false;
    break;
  case 0:
    {
      PermGroup::ConstructionMethod constr;
      if (schreier_sims == "standard")
        constr = PermGroup::SCHREIER_SIMS;
      else if (schreier_sims == "random")
        constr = PermGroup::SCHREIER_SIMS_RANDOM;
      else {
        error("unreachable");
        return false;
      }

      PermGroup::TransversalStorageMethod transv;
      if (transversal_storage == "explicit")
        transv = PermGroup::EXPLICIT_TRANSVERSALS;
      else if (transversal_storage == "schreier-trees")
        transv = PermGroup::SCHREIER_TREES;
      else if (transversal_storage == "shallow-schreier-trees")
        transv = PermGroup::SHALLOW_SCHREIER_TREES;
      else {
        error("unreachable");
        return false;
      }

      PermGroup g(degree, gens, constr, transv);
    }
    break;
  default:
    if (!time_child(child, t, &t_acc))
      return false;
  }

  return true;
}

bool run_gap(std::string const &generators, double *t)
{
  static double t_acc = 0.0;

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

  f << "StabChain(Group(" + generators + "));\n";
  f.flush();

  pid_t child;
  switch((child = fork())) {
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
    }
    break;
  default:
    if (!time_child(child, t, &t_acc))
      return false;
  }

  return true;
}

// statistics

void get_statistics(std::vector<double> const &ts, double *mean, double *stddev) {
  *mean = std::accumulate(ts.begin(), ts.end(), 0.0) / ts.size();

  std::vector<double> d(ts.size());
  std::transform(ts.begin(), ts.end(), d.begin(),
                 [mean](double t) { return t - *mean; });

  *stddev = std::sqrt(
    std::inner_product(d.begin(), d.end(), d.begin(), 0.0) / ts.size());
}

}

int main(int argc, char **argv)
{
  progname = basename(argv[0]);

  struct option long_options[] = {
    {"help",                  no_argument,       0,       'h'},
    {"schreier-sims",         required_argument, 0,       's'},
    {"transversal-storage",   required_argument, 0,       't'},
    {"num-runs",              required_argument, 0,       'n'},
    {"verbose",               no_argument,       0,       'v'},
    {nullptr,                 0,                 nullptr,  0 }
  };

  std::string schreier_sims;
  std::string transversal_storage;
  std::string groups;
  int num_runs = 1;
  bool verbose = false;

  // TODO: detect superfluous non-argument options

  int c;
  while ((c = getopt_long(argc, argv, "hs:t:n:v", long_options, nullptr)) != -1) {
    switch(c) {
    case 'h':
      usage(std::cout);
      return EXIT_SUCCESS;
    case 's':
      {
        if (!valid_choice(schreier_sims_choices, optarg)) {
          error("invalid schreier-sims variant specified");
          return EXIT_FAILURE;
        }
        schreier_sims = optarg;
      }
      break;
    case 't':
      {
        if (!valid_choice(transversal_storage_choices, optarg)) {
          error("invalid transversal-storage variant specified");
          return EXIT_FAILURE;
        }
        transversal_storage = optarg;
      }
      break;
    case 'n':
      {
        bool failed = false;

        std::size_t idx;
        try {
          num_runs = std::stoi(optarg, &idx);
          failed = idx != strlen(optarg);
        } catch (...) {
          failed = true;
        }

        if (failed) {
          error("invalid argument to --num-runs");
          return EXIT_FAILURE;
        }
      }
      break;
    case 'v':
      verbose = true;
      break;
    default:
      break;
    }
  }

  if (optind == argc) {
    usage(std::cerr);
    error("GROUPS argument is mandatory");
    return EXIT_FAILURE;
  }

  groups = argv[optind];

  if (schreier_sims.empty()) {
    usage(std::cerr);
    error("--schreier-sims option is mandatory");
    return EXIT_FAILURE;
  }

  if (schreier_sims != "gap" && transversal_storage.empty()) {
    usage(std::cerr);
    error("--transversal-storage option is mandatory when not using gap");
    return EXIT_FAILURE;
  }

  std::ifstream f(groups);
  if (f.fail()) {
    error("failed to open " + groups);
    return EXIT_FAILURE;
  }

  std::regex re("degree:(\\d+),order:(\\d+),gens:(.*)");
  std::smatch m;

  std::string line;
  int lineno = 1;
  while (std::getline(f, line)) {
    if (!std::regex_match(line, m, re)) {
      error("failed to parse line no.", lineno, "in", groups);
      return EXIT_FAILURE;
    }

    if (verbose) {
      std::cerr << "profiling group " << lineno
                << " with degree " << m[1]
                << " and order " << m[2]
                << '\n';
    }

    std::vector<double> ts;
    for (int r = 0; r < num_runs; ++r) {
      if (verbose)
        std::cerr << "run " << r + 1 << '/' << num_runs << '\n';

      double t;
      if (schreier_sims == "gap") {
        if (!run_gap(m[3], &t)) {
          error("failed to run gap");
          return EXIT_FAILURE;
        }
      } else {
        if (!run_cpp(m[3], schreier_sims, transversal_storage, &t))
          return EXIT_FAILURE;
      }

      ts.push_back(t);
    }

    double t_mean, t_stddev;
    get_statistics(ts, &t_mean, &t_stddev);

    std::cout.precision(3);
    std::cout << "mean: " << std::fixed << t_mean << "s, "
              << "stddev: " << std::fixed << t_stddev << "s\n";

    ++lineno;
  }

  if (f.bad()) {
    error("failed to read " + groups);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
