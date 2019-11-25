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

#include "arch_graph_system.h"
#include "perm_set.h"
#include "profile_args.h"
#include "profile_parse.h"
#include "profile_run.h"
#include "profile_timer.h"
#include "profile_utility.h"
#include "task_mapping.h"
#include "task_orbits.h"

namespace
{

std::string progname;

void usage(std::ostream &s)
{
  char const *opts[] = {
    "[-h|--help]",
    "-i|--implementation {gap|mpsym}",
    "[-a|--approximate]",
    "[--realtime-clock]",
    "[-v|--verbose]",
    "GROUP",
    "TASK_ALLOCATIONS"
  };

  s << "usage: " << progname << '\n';
  for (char const *opt : opts)
    s << "  " << opt << '\n';
}

std::string map_tasks_gap(unsigned degree,
                          std::string const &generators,
                          std::string const &task_allocations,
                          bool verbose)
{
  (void) degree;

  std::stringstream ss;

  ss << "LoadPackage(\"orb\");\n";

  ss << "automorphisms:=Group(" << generators << ");\n";

  ss << "task_allocations:=[\n";
  ss << task_allocations;
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

  return ss.str();
}

void map_tasks_mpsym(bool approximate,
                     unsigned degree,
                     cgtl::PermSet const &generators,
                     std::vector<cgtl::TaskAllocation> const &task_allocations,
                     bool verbose)
{
  using cgtl::ArchGraphSystem;
  using cgtl::PermGroup;
  using cgtl::TaskOrbits;

  ArchGraphSystem ag(PermGroup(degree, generators));

  TaskOrbits task_orbits;
  for (auto const &ta : task_allocations)
    task_orbits.insert(ag.mapping({ta, 0u, approximate}));

  if (verbose)
    std::cout << "Found " << task_orbits.num_orbits() << " equivalence classes"
              << std::endl;
}

double run(VariantOption const &library_impl,
           bool approximate,
           unsigned degree,
           std::string const &generators_str,
           std::string const &task_allocations_str,
           bool verbose)
{
  double t;
  if (library_impl.is("gap")) {
    run_gap(
      map_tasks_gap(degree,
                    parse_generators_gap(generators_str),
                    parse_task_allocations_gap(task_allocations_str),
                    verbose),
      &t);

  } else if (library_impl.is("mpsym")) {
    run_cpp([&]{
      map_tasks_mpsym(approximate,
                      degree,
                      parse_generators_mpsym(generators_str),
                      parse_task_allocations_mpsym(task_allocations_str),
                      verbose);
      }, &t);
  }

  return t;
}

void profile(VariantOption const &library_impl,
             bool approximate,
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

  double t = run(library_impl,
                 approximate,
                 degree,
                 generators_str,
                 task_allocations_str,
                 verbose);

  std::cout << "Runtime: " << std::scientific << t << std::endl;
}

} // namespace

int main(int argc, char **argv)
{
  progname = basename(argv[0]);

  struct option long_options[] = {
    {"help",                no_argument,       0,       'h'},
    {"implementation",      required_argument, 0,       'i'},
    {"approximate",         no_argument,       0,       'a'},
// TODO: check correctness
    {"realtime-clock",      no_argument,       0,        1 },
    {"verbose",             no_argument,       0,       'v'},
    {nullptr,               0,                 nullptr,  0 }
  };

  VariantOption library_impl({"gap", "mpsym"});

  bool approximate = false;

  bool verbose = false;

  for (;;) {
    int c = getopt_long(argc, argv, "hi:av", long_options, nullptr);
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
      case 'a':
        approximate = true;
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

  CHECK_OPTION((!approximate || !library_impl.is("gap")),
               "--approximate not supported when using gap");

  try {
    profile(library_impl,
            approximate,
            group_stream,
            task_allocations_stream,
            verbose);
  } catch (std::exception const &e) {
    error("profiling failed", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
