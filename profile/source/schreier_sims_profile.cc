#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>

#include <getopt.h>
#include <libgen.h>


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

void run_cpp(std::string const &generators,
             std::string const &schreier_sims,
             std::string const &transversal_storage)
{
  throw std::logic_error("TODO");
}

void run_gap(std::string const &generators)
{
   // TODO
   std::cout << "StabChain(Group(" << generators << "))\n";
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
    {nullptr,                 0,                 nullptr,  0 }
  };

  std::string schreier_sims;
  std::string transversal_storage;
  std::string groups;
  int num_runs;

  while (optind < argc) {
    int option_index = 0;
    int c = getopt_long(argc, argv, "hs:t:n:", long_options, &option_index);

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
    default:
      if (!groups.empty()) {
        usage(std::cerr);
        return EXIT_FAILURE;
      }
      groups = argv[optind++];
      break;
    }
  }

  if (schreier_sims.empty()) {
    usage(std::cerr);
    error("--schreier-sims option is mandatory");
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

    if (schreier_sims == "gap")
      run_gap(m[3]);
    else
      throw std::logic_error("TODO");

    ++lineno;
  }

  if (f.bad()) {
    error("failed to read " + groups);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
