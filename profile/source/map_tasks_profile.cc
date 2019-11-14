#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <vector>

#include <getopt.h>
#include <libgen.h>

#include "profile_args.h"
#include "profile_parse.h"
#include "profile_run.h"
#include "profile_timer.h"
#include "profile_utility.h"

namespace
{

std::string progname;

void usage(std::ostream &s)
{
  char const *opts[] = {
    "[-h|--help]",
    "-i|--implementation {gap|mpsym}",
    "GROUP",
    "TASK_ALLOCATIONS"
  };

  s << "usage: " << progname << '\n';
  for (char const *opt : opts)
    s << "  " << opt << '\n';
}

void run(VariantOption const &library_impl,
         unsigned degree,
         std::string const &generators_str,
         std::string const &task_allocations_str,
         bool verbose)
{
  (void)degree;

  double t;
  if (library_impl.is("gap")) {
    auto generators(parse_generators_gap(generators_str));

    std::stringstream ss;

    ss << "LoadPackage(\"orb\");\n";

    ss << "automorphisms:=Group(" << parse_generators_gap(generators_str) << ");\n";

    ss << "task_allocations:=[\n";
    ss << parse_task_allocations_gap(task_allocations_str);
    ss << "];\n";

    ss << "orbit_representatives:=HTCreate([1,2,3]);\n";
    ss << "orbit_options:=rec(lookingfor:=orbit_representatives);\n";

    ss << "for task_allocation in task_allocations do\n";
    ss << "  orbit:=Orb(automorphisms, task_allocation, OnTuples, orbit_options);\n";
    ss << "  if not PositionOfFound(orbit) then\n";
    ss << "    orbit_representative:=Elements(Enumerate(orbit))[1];\n";
    ss << "    HTAdd(orbit_representatives, orbit_representative, true);\n";
    ss << "  fi;\n";
    ss << "od;\n";

    if (verbose)
      ss << "Print(\"Found \", orbit_representatives!.nr, \" equivalence classes\\n\");\n";

    run_gap(ss.str(), &t);

  } else if (library_impl.is("mpsym")) {
    throw std::logic_error("TODO");
  }

  std::cout << "Runtime: " << std::scientific << t << std::endl;
}

void profile(VariantOption const &library_impl,
             std::ifstream &group_stream,
             std::ifstream &task_allocations_stream,
             bool verbose)
{
  std::string line;
  if (!std::getline(group_stream, line))
    throw std::runtime_error("failed to read group file");

  unsigned degree;
  unsigned order;
  std::string generators_str;

  std::tie(degree, order, generators_str) = parse_group(line);

  if (verbose) {
    std::cout << ">>> Using automorphism group "
              << " with degree " << degree
              << ", order " << order
              << " and generators " << generators_str
              << std::endl;
  }

  std::string task_allocations_str(
    (std::istreambuf_iterator<char>(task_allocations_stream)),
    std::istreambuf_iterator<char>());

  run(library_impl,
      degree,
      generators_str,
      task_allocations_str,
      verbose);
}

} // namespace

int main(int argc, char **argv)
{
  progname = basename(argv[0]);

  struct option long_options[] = {
    {"help",                no_argument,       0,       'h'},
    {"implementation",      required_argument, 0,       'i'},
    {"realtime-clock",      no_argument,       0,        1 },
    {"verbose",             no_argument,       0,       'v'},
    {nullptr,               0,                 nullptr,  0 }
  };

  VariantOption library_impl({"gap", "mpsym"});

  bool verbose = false;

  for (;;) {
    int c = getopt_long(argc, argv, "hi:v", long_options, nullptr);
    if (c == -1)
      break;

    try {
      switch(c) {
      case 'h':
        usage(std::cout);
        return EXIT_SUCCESS;
      case 'i':
        library_impl.set(optarg);
        break;
      case 'v':
        verbose = true;
        break;
      case 1:
        timer_enable_realtime();
        break;
      default:
        return EXIT_FAILURE;
      }
    } catch (std::invalid_argument const &) {
      error("invalid option argument");
      return EXIT_FAILURE;
    }
  }

  CHECK_OPTION(library_impl.is_set(),
               "--implementation option is mandatory");

  CHECK_FILE_ARGUMENT("GROUP");

  std::ifstream group_stream(argv[optind++]);

  CHECK_FILE_ARGUMENT("TASK_ALLOCATIONS");

  std::ifstream task_allocations_stream(argv[optind]);

  try {
    profile(library_impl, group_stream, task_allocations_stream, verbose);
  } catch (std::exception const &e) {
    error("profiling failed", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
