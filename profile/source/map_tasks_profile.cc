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
#include "dump.h"
#include "perm_set.h"
#include "task_mapping.h"
#include "task_orbits.h"
#include "timer.h"

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

struct ProfileOptions
{
  VariantOption library;
  bool approximate;
  bool verbose;
};

std::string map_tasks_gap(std::string const &generators,
                          std::string const &task_allocations,
                          ProfileOptions const &options)
{
  std::stringstream ss;

  ss << "LoadPackage(\"orb\");\n";

  ss << "automorphisms:=Group(" << generators << ");\n";

  ss << "task_allocations:=[\n";
  ss << task_allocations;
  ss << "];\n";

  ss << "orbit_representatives:=[];\n";
  ss << "orbit_representatives_hash:=HTCreate([1,2,3]);\n";

  ss << "orbit_options:=rec(lookingfor:=orbit_representatives_hash);\n";

  ss << "n:=1;\n";
  ss << "for task_allocation in task_allocations do\n";

  if (options.verbose)
    ss << "  Print(\"INFO: Mapping task \", n, \" of \", "
          "        Length(task_allocations), \"\\r\");\n";

  ss << "  orbit:=Orb(automorphisms, task_allocation, OnTuples, orbit_options);\n";

  ss << "  if not PositionOfFound(orbit) then\n";
  ss << "    orbit_repr:=Elements(Enumerate(orbit))[1];\n";

  ss << "    if HTAdd(orbit_representatives_hash, orbit_repr, true) <> fail then\n";
  ss << "      Append(orbit_representatives, [orbit_repr]);\n";
  ss << "    fi;\n";
  ss << "  fi;\n";

  ss << "  n:=n+1;\n";
  ss << "od;\n";

  if (options.verbose) {
    ss << "Print(\"\\nINFO: Found \", Length(orbit_representatives), "
          "      \" equivalence classes:\\n\");\n";
    ss << "for orbit_repr in orbit_representatives do\n";
    ss << "  Print(\"INFO: \", orbit_repr, \"\\n\");\n";
    ss << "od;\n";
  }

  return ss.str();
}

void map_tasks_mpsym(cgtl::PermSet const &generators,
                     std::vector<cgtl::TaskAllocation> const &task_allocations,
                     ProfileOptions const &options)
{
  using cgtl::ArchGraphSystem;
  using cgtl::PermGroup;
  using cgtl::TaskOrbits;

  ArchGraphSystem ag(PermGroup(generators.degree(), generators));

  TaskOrbits task_orbits;
  for (auto i = 0u; i < task_allocations.size(); ++i) {
    if (options.verbose)
      progress("Mapping task", i + 1u, "of", task_allocations.size());

    task_orbits.insert(ag.mapping({task_allocations[i], 0u, options.approximate}));
  }

  if (options.verbose) {
    progress_done();

    info("Found", task_orbits.num_orbits(), "equivalence classes:");

    for (auto repr : task_orbits)
      info(dump::dump(repr));

    info("Timer dumps:");
    Timer_dump(options.approximate ? "map approx" : "map bruteforce");
  }
}

double run(std::string const &generators_str,
           std::string const &task_allocations_str,
           ProfileOptions const &options)
{
  double t;
  if (options.library.is("gap")) {
    run_gap(
      map_tasks_gap(parse_generators_gap(generators_str),
                    parse_task_allocations_gap(task_allocations_str),
                    options),
      &t);

  } else if (options.library.is("mpsym")) {
    run_cpp([&]{
      map_tasks_mpsym(parse_generators_mpsym(generators_str),
                      parse_task_allocations_mpsym(task_allocations_str),
                      options);
      }, &t);
  }

  return t;
}

void profile(std::ifstream &group_stream,
             std::ifstream &task_allocations_stream,
             ProfileOptions const &options)
{
  std::string line;
  if (!std::getline(group_stream, line))
    throw std::runtime_error("failed to read group file");

  unsigned degree;
  unsigned order;
  std::string generators_str;

  std::tie(degree, order, generators_str) = parse_group(line);

  if (options.verbose) {
    info("Using automorphism group with degree", degree,
         "order", order,
         "and generators", generators_str);
  }

  std::string task_allocations_str(
    (std::istreambuf_iterator<char>(task_allocations_stream)),
    std::istreambuf_iterator<char>());

  double t = run(generators_str,
                 task_allocations_str,
                 options);

  result("Runtime:", t);
}

} // namespace

int main(int argc, char **argv)
{
  progname = basename(argv[0]);

  struct option long_options[] = {
    {"help",                no_argument,       0,       'h'},
    {"implementation",      required_argument, 0,       'i'},
    {"approximate",         no_argument,       0,       'a'},
    {"realtime-clock",      no_argument,       0,        1 },
    {"verbose",             no_argument,       0,       'v'},
    {nullptr,               0,                 nullptr,  0 }
  };

  VariantOption library({"gap", "mpsym"});

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
        library.set(optarg);
        break;
      case 'a':
        approximate = true;
        break;
      case 'v':
        verbose = true;
        Timer::enabled = true;
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

  CHECK_OPTION(library.is_set(),
               "--implementation option is mandatory");

  CHECK_FILE_ARGUMENT(group_stream, "GROUP");

  CHECK_FILE_ARGUMENT(task_allocations_stream, "TASK_ALLOCATIONS");

  CHECK_OPTION((!approximate || !library.is("gap")),
               "--approximate not supported when using gap");

  try {
    profile(group_stream,
            task_allocations_stream,
            {library, approximate, verbose});
  } catch (std::exception const &e) {
    error("profiling failed", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
