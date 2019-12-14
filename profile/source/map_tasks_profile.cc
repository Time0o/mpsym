#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
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
#include "profile_generate.h"
#include "profile_parse.h"
#include "profile_read.h"
#include "profile_run.h"
#include "profile_timer.h"
#include "profile_util.h"

namespace
{

std::string progname;

void usage(std::ostream &s)
{
  char const *opts[] = {
    "[-h|--help]",
    "-i|--implementation {gap|mpsym|mpsym_approx}",
    "-g|--groups GROUPS",
    "[-t|--task-allocations TASK_ALLOCATIONS]",
    "[--num-tasks NUM_TASKS]",
    "[--num-task-allocations NUM_TASK_ALLOCATIONS]",
    "[--realtime-clock]",
    "[-v|--verbose]"
  };

  s << "usage: " << progname << '\n';
  for (char const *opt : opts)
    s << "  " << opt << '\n';
}

struct ProfileOptions
{
  VariantOption library{"gap", "mpsym", "mpsym_approx"};
  unsigned num_tasks = 0u;
  unsigned num_task_allocations = 0u;
  bool verbose = false;
};

std::string map_tasks_gap(gap::PermSet const &generators,
                          gap::TaskAllocationVector const &task_allocation_vector,
                          ProfileOptions const &options)
{
  std::stringstream ss;

  ss << "LoadPackage(\"orb\");\n";

  ss << "automorphisms:=Group(" << generators.permutations << ");\n";

  ss << "task_allocations:=[\n";
  ss << task_allocation_vector.task_allocations;
  ss << "];\n";

  ss << "orbit_representatives:=[];\n";
  ss << "orbit_representatives_hash:=HTCreate([1,2,3]);\n";

  ss << "orbit_options:=rec(lookingfor:=orbit_representatives_hash);\n";

  ss << "n:=1;\n";
  ss << "for task_allocation in task_allocations do\n";

  if (options.verbose)
    ss << "  Print(\"INFO: Mapping task \", n, \" of \", "
          "        Length(task_allocations), \"\\r\\c\");\n";

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
                     cgtl::TaskAllocationVector const &task_allocation_vector,
                     ProfileOptions const &options)
{
  using cgtl::ArchGraphSystem;
  using cgtl::PermGroup;
  using cgtl::TaskOrbits;

  ArchGraphSystem ag(PermGroup(generators.degree(), generators));
  auto task_allocations(task_allocation_vector.task_allocations);

  TaskOrbits task_orbits;
  for (auto i = 0u; i < task_allocations.size(); ++i) {
    if (options.verbose)
      progress("Mapping task", i + 1u, "of", task_allocations.size());

    auto mapping(ag.mapping(
      {task_allocations[i], 0u, options.library.is("mpsym_approx")}));

    task_orbits.insert(mapping);
  }

  if (options.verbose) {
    progress_done();

    info("Found", task_orbits.num_orbits(), "equivalence classes:");

    for (auto const &repr : task_orbits)
      info(dump::dump(repr));

    info("Timer dumps:");
    Timer_dump(options.library.is("mpsym_approx") ? "map approx"
                                                  : "map bruteforce");
  }
}

double run(std::string const &generators,
           std::string const &task_allocations,
           ProfileOptions const &options)
{
  double t;
  if (options.library.is("gap")) {
    auto generators_gap(parse_generators_gap(generators));
    auto task_allocations_gap(parse_task_allocations_gap(task_allocations));

    if (task_allocations_gap.max_pe > generators_gap.degree)
      throw std::invalid_argument("pe index out of range");

    run_gap(map_tasks_gap(generators_gap, task_allocations_gap, options), &t);

  } else if (options.library.is("mpsym") || options.library.is("mpsym_approx")) {
    auto generators_mpsym(parse_generators_mpsym(generators));
    auto task_allocations_mpsym(parse_task_allocations_mpsym(task_allocations));

    if (task_allocations_mpsym.max_pe > generators_mpsym.degree())
      throw std::invalid_argument("pe index out of range");

    run_cpp([&]{
      map_tasks_mpsym(generators_mpsym, task_allocations_mpsym, options);
    }, &t);
  }

  return t;
}

void profile(Stream &groups_stream,
             Stream &task_allocations_stream,
             ProfileOptions const &options)
{
  if (options.verbose)
    info("Implementation:", options.library.get());

  std::string task_allocations;

  if (task_allocations_stream.valid)
    task_allocations = read_file(task_allocations_stream.stream);

  foreach_line(groups_stream.stream, [&](std::string const &line, unsigned lineno){
    auto group(parse_group(line));

    if (!task_allocations_stream.valid)
      task_allocations = generate_task_allocations(
        group.degree, options.num_tasks, options.num_task_allocations);

    if (options.verbose) {
      info("Using automorphism group", lineno,
           "with degree", group.degree,
           "and generators", group.generators);
    }

    double t = run(group.generators, task_allocations, options);

    result("Runtime:", t);
  });
}

} // namespace

int main(int argc, char **argv)
{
  progname = basename(argv[0]);

  struct option long_options[] = {
    {"help",                 no_argument,       0,       'h'},
    {"implementation",       required_argument, 0,       'i'},
    {"groups",               required_argument, 0,       'g'},
    {"task-allocations",     required_argument, 0,       't'},
    {"num-tasks",            required_argument, 0,        2 },
    {"num-task-allocations", required_argument, 0,        3 },
    {"realtime-clock",       no_argument,       0,        4 },
    {"verbose",              no_argument,       0,       'v'},
    {nullptr,                0,                 nullptr,  0 }
  };

  ProfileOptions options;

  Stream groups_stream;
  Stream task_allocations_stream;

  for (;;) {
    int c = getopt_long(argc, argv, "hi:g:t:v", long_options, nullptr);
    if (c == -1)
      break;

    try {
      switch(c) {
      case 'h':
        usage(std::cout);
        return EXIT_SUCCESS;
      case 'i':
        options.library.set(optarg);
        break;
      case 'g':
        OPEN_STREAM(groups_stream, optarg);
        break;
      case 't':
        OPEN_STREAM(task_allocations_stream, optarg);
        break;
      case 2:
        options.num_tasks = stox<unsigned>(optarg);
        break;
      case 3:
        options.num_task_allocations = stox<unsigned>(optarg);
        break;
      case 'v':
        options.verbose = true;
        Timer::enabled = true;
        break;
      case 4:
        timer_realtime_enable();
        break;
      default:
        return EXIT_FAILURE;
      }
    } catch (std::invalid_argument const &) {
      error("invalid option argument");
      return EXIT_FAILURE;
    }
  }

  CHECK_OPTION(options.library.is_set(), "--implementation option is mandatory");

  CHECK_OPTION(groups_stream.valid, "--groups option is mandatory");

  if (task_allocations_stream.valid) {
    if (options.num_tasks > 0 || options.num_task_allocations > 0)
       warning("task allocations explicitly given,"
               " --num-tasks, --num-task-allocations ignored");

  } else {
    CHECK_OPTION(options.num_tasks > 0 || options.num_task_allocations > 0,
                 "task allocations not explicitly given,"
                 " --num-tasks, --num-task-allocations missing");
  }

  try {
    profile(groups_stream, task_allocations_stream, options);
  } catch (std::exception const &e) {
    error("profiling failed:", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
