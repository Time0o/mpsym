#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <getopt.h>
#include <libgen.h>

#include "arch_graph_system.h"
#include "dump.h"
#include "perm_set.h"
#include "task_allocation.h"
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
    "-i|--implementation {gap|mpsym}",
    "-m|--mapping-method {iterate|local_search|orbits}",
    "[--mapping-options {dont_match_reprs}]",
    "-g|--groups GROUPS",
    "[-t|--task-allocations TASK_ALLOCATIONS]",
    "[--num-tasks NUM_TASKS]",
    "[--num-task-allocations NUM_TASK_ALLOCATIONS]",
    "[--check-accuracy]",
    "[--realtime-clock]",
    "[-v|--verbose]",
    "[--show-gap-errors]"
  };

  s << "usage: " << progname << '\n';
  for (char const *opt : opts)
    s << "  " << opt << '\n';
}

struct ProfileOptions
{
  VariantOption library{"gap", "mpsym"};
  VariantOption mapping_method{"iterate", "local_search", "orbits"};
  VariantOptionSet mapping_options{"dont_match_reprs"};
  unsigned num_tasks = 0u;
  unsigned num_task_allocations = 0u;
  bool check_accuracy = false;
  int verbosity = 0;
  bool show_gap_errors = false;
};

std::string map_tasks_gap_iterate(ProfileOptions const &options)
{
  std::stringstream ss;

  ss << "orbit_repr:=task_allocation;\n";
  ss << "orbit_repr_new:=true;\n";

  ss << "for element in automorphisms do\n";
  ss << "  permuted:=OnTuples(task_allocation, element);\n";

  if (options.mapping_options.is_set("dont_match_reprs")) {
    ss << "  if permuted < orbit_repr then\n";
    ss << "    orbit_repr:=permuted;\n";
    ss << "  fi;\n";
    ss << "od;\n";

    ss << "if HTAdd(orbit_representatives_hash, orbit_repr, true) <> fail then\n";
    ss << "  Append(orbit_representatives, [orbit_repr]);\n";
    ss << "fi;\n";

  } else {
    ss << "  if HTValue(orbit_representatives_hash, permuted) <> fail then\n";
    ss << "    orbit_repr_new:=false;\n";
    ss << "    break;\n";
    ss << "  elif permuted < orbit_repr then\n";
    ss << "    orbit_repr:=permuted;\n";
    ss << "  fi;\n";
    ss << "od;\n";

    ss << "if orbit_repr_new then\n";
    ss << "  HTAdd(orbit_representatives_hash, orbit_repr, true);\n";
    ss << "  Append(orbit_representatives, [orbit_repr]);\n";
    ss << "fi;\n";

  }

  return ss.str();
}

std::string map_tasks_gap_orbits(ProfileOptions const &options)
{
  std::stringstream ss;

  if (options.mapping_options.is_set("dont_match_reprs")) {
    ss << "orbit:=Orb(automorphisms, task_allocation, OnTuples);\n";
    ss << "orbit_repr:=Elements(Enumerate(orbit))[1];\n";

    ss << "if HTAdd(orbit_representatives_hash, orbit_repr, true) <> fail then\n";
    ss << "  Append(orbit_representatives, [orbit_repr]);\n";
    ss << "fi;\n";

  } else {
    ss << "orbit_options:=rec(lookingfor:=orbit_representatives_hash);\n";
    ss << "orbit:=Orb(automorphisms, task_allocation, OnTuples, orbit_options);\n";

    ss << "if not PositionOfFound(orbit) then\n";
    ss << "  orbit_repr:=Elements(Enumerate(orbit))[1];\n";

    ss << "  HTAdd(orbit_representatives_hash, orbit_repr, true);\n";
    ss << "  Append(orbit_representatives, [orbit_repr]);\n";
    ss << "fi;\n";
  }

  return ss.str();
}

std::string map_tasks_gap(
  gap::PermSet const &generators,
  gap::TaskAllocationVector const &task_allocation_vector,
  ProfileOptions const &options)
{
  std::stringstream ss;

  // load the "orb" package containing orbit enumeration and hashing functions
  ss << "LoadPackage(\"orb\");\n";

  // construct the automorphism group
  ss << "automorphisms:=Group(" << generators.permutations << ");\n";

  // construct the vector of task allocations to be mapped
  ss << "task_allocations:=[\n";
  ss << task_allocation_vector.task_allocations;
  ss << "];\n";

  ss << "orbit_representatives:=[];\n";
  ss << "orbit_representatives_hash:=HTCreate([1,2,3]);\n";

  // map tasks allocations one by one
  ss << "n:=1;\n";
  ss << "for task_allocation in task_allocations do\n";

  // display progress
  if (options.verbosity > 0) {
    ss << "  Print(\"DEBUG: Mapping task \", n, \" of \", "
          "        Length(task_allocations), \"\\r\\c\");\n";
  }

  // concrete mapping code depending on chosen implementation
  if (options.mapping_method.is("iterate"))
    ss << map_tasks_gap_iterate(options);
  else if (options.mapping_method.is("orbits"))
    ss << map_tasks_gap_orbits(options);
  else
    throw std::logic_error("unreachable");

  ss << "  n:=n+1;\n";
  ss << "od;\n";

  // display orbit representatives found
  if (options.check_accuracy || options.verbosity > 0) {
    ss << "Print(\"\\nDEBUG: => Found \", Length(orbit_representatives), "
          "      \" orbit representatives\\n\");\n";

    if (options.check_accuracy || options.verbosity > 1) {
      ss << "for orbit_repr in orbit_representatives do\n";
      ss << "  Print(\"DEBUG: \", orbit_repr, \"\\n\");\n";
      ss << "od;\n";
    }
  }

  return ss.str();
}

cgtl::TaskOrbits map_tasks_mpsym(
  cgtl::PermSet const &generators,
  cgtl::TaskAllocationVector const &task_allocation_vector,
  ProfileOptions const &options)
{
  using cgtl::ArchGraphSystem;
  using cgtl::PermGroup;
  using cgtl::TaskOrbits;

  ArchGraphSystem ag(PermGroup(generators.degree(), generators));

  auto task_allocations(task_allocation_vector.task_allocations);

  // setup mapping options
  ArchGraphSystem::MappingOptions mapping_options;

  if (options.mapping_method.is("iterate"))
    mapping_options.method = ArchGraphSystem::MappingMethod::ITERATE;
  else if (options.mapping_method.is("local_search"))
    mapping_options.method = ArchGraphSystem::MappingMethod::LOCAL_SEARCH;
  else if (options.mapping_method.is("orbits"))
    mapping_options.method = ArchGraphSystem::MappingMethod::ORBITS;
  else
    throw std::logic_error("unreachable");

  if (options.mapping_options.is_set("dont_match_reprs"))
      mapping_options.match_reprs = false;

  // perform mappings
  TaskOrbits task_orbits;
  for (auto i = 0u; i < task_allocations.size(); ++i) {
    if (options.verbosity > 0)
      debug_progress("Mapping task", i + 1u, "of", task_allocations.size());

    auto mapping(ag.mapping(
      task_allocations[i], 0u, &mapping_options, &task_orbits));
  }

  if (options.verbosity > 0) {
    debug_progress_done();

    debug("=> Found", task_orbits.num_orbits(), "orbit representatives");
    if (options.verbosity > 1) {
      for (auto const &repr : task_orbits)
        debug(DUMP(repr));
    }

    debug("Timer dumps:");
    if (options.mapping_method.is("iterate"))
      debug_timer_dump("map bruteforce iterate");
    else if (options.mapping_method.is("local_search"))
      debug_timer_dump("map approx local search");
    else if (options.mapping_method.is("orbits"))
      debug_timer_dump("map bruteforce orbits");
  }

  return task_orbits;
}

void map_tasks_gap_wrapper(std::string const &generators,
                           std::string const &task_allocations,
                           ProfileOptions const &options,
                           cgtl::TaskOrbits *task_orbits,
                           double *t)
{
  using cgtl::TaskOrbits;

  // parse input

  auto generators_gap(parse_generators_gap(generators));
  auto task_allocations_gap(parse_task_allocations_gap(task_allocations));

  if (task_allocations_gap.max_pe > generators_gap.degree)
    throw std::invalid_argument("pe index out of range");

  // run gap script

  auto gap_script(map_tasks_gap(generators_gap,
                                task_allocations_gap,
                                options));

  auto gap_output(run_gap(gap_script,
                          options.check_accuracy,
                          !options.show_gap_errors,
                          t));

  // parse output

  if (options.check_accuracy) {
    auto representatives_gap(parse_task_allocations_gap_to_mpsym(gap_output));

    task_orbits->insert_all(representatives_gap.task_allocations.begin(),
                            representatives_gap.task_allocations.end());
  }
}

void map_tasks_mpsym_wrapper(std::string const &generators,
                             std::string const &task_allocations,
                             ProfileOptions const &options,
                             cgtl::TaskOrbits *task_orbits,
                             double *t)
{
  using cgtl::TaskOrbits;

  // parse input

  auto generators_mpsym(parse_generators_mpsym(generators));
  auto task_allocations_mpsym(parse_task_allocations_mpsym(task_allocations));

  if (task_allocations_mpsym.max_pe > generators_mpsym.degree())
    throw std::invalid_argument("pe index out of range");

  // run task mapping code

  *task_orbits = run_cpp([&]{
    return map_tasks_mpsym(generators_mpsym,
                           task_allocations_mpsym,
                           options);
  }, t);
}

void check_accuracy(cgtl::TaskOrbits const &task_orbits_mpsym,
                    cgtl::TaskOrbits const &task_orbits_gap,
                    ProfileOptions const &options)
{
  using cgtl::TaskAllocation;

  if (task_orbits_mpsym == task_orbits_gap) {
    info("Orbit representatives match");
    return;
  }

  info("Orbit representatives do not match:");

  // construct representative sets
  std::set<TaskAllocation> reprs_mpsym, reprs_gap, reprs_missing, reprs_extra;

  for (auto const &repr : task_orbits_mpsym)
    reprs_mpsym.insert(repr);

  for (auto const &repr : task_orbits_gap)
    reprs_gap.insert(repr);

  // find missing/extra representatives
  if (!reprs_gap.empty()) {
    for (auto const &repr : reprs_gap) {
      if (reprs_mpsym.find(repr) == reprs_mpsym.end())
        reprs_missing.insert(repr);
    }

    info("=>", reprs_missing.size(), "Missing orbit representatives");
    if (options.verbosity > 1) {
      for (auto const &repr : reprs_missing)
        info(DUMP(repr));
    }
  }

  if (!reprs_mpsym.empty()) {
    for (auto const &repr : reprs_mpsym) {
      if (reprs_gap.find(repr) == reprs_gap.end())
        reprs_extra.insert(repr);
    }

    info("=>", reprs_extra.size(), "Missing orbit representatives");
    if (options.verbosity > 1) {
      for (auto const &repr : reprs_extra)
        info(DUMP(repr));
    }
  }
}

double run(std::string const &generators,
           std::string const &task_allocations,
           ProfileOptions const &options)
{
  using cgtl::TaskOrbits;

  double t;

  if (options.library.is("gap")) {
    map_tasks_gap_wrapper(generators, task_allocations, options, nullptr, &t);

  } else if (options.library.is("mpsym")) {
    TaskOrbits task_orbits_mpsym, task_orbits_gap;

    map_tasks_mpsym_wrapper(generators,
                            task_allocations,
                            options,
                            &task_orbits_mpsym,
                            &t);

    if (options.check_accuracy) {
      info("Checking accuracy...");

      map_tasks_gap_wrapper(generators,
                            task_allocations,
                            options,
                            &task_orbits_gap,
                            nullptr);

      check_accuracy(task_orbits_mpsym, task_orbits_gap, options);
    }
  }

  return t;
}

void profile(Stream &groups_stream,
             Stream &task_allocations_stream,
             ProfileOptions const &options)
{
  if (options.verbosity > 0)
    debug("Implementation:", options.library.get());

  std::string task_allocations;

  if (task_allocations_stream.valid)
    task_allocations = read_file(task_allocations_stream.stream);

  foreach_line(groups_stream.stream, [&](std::string const &line, unsigned lineno){
    auto group(parse_group(line));

    if (!task_allocations_stream.valid)
      task_allocations = generate_task_allocations(
        group.degree, options.num_tasks, options.num_task_allocations);

    if (options.verbosity > 0) {
      info("Using automorphism group", lineno,
           "with degree", group.degree,
           "and generators", group.generators);
    } else {
      info("Using automorphism group", lineno,
           "with degree", group.degree);
    }

    double t = run(group.generators, task_allocations, options);

    result("Runtime:", t, "s");
  });
}

} // namespace

int main(int argc, char **argv)
{
  progname = basename(argv[0]);

  struct option long_options[] = {
    {"help",                 no_argument,       0,       'h'},
    {"implementation",       required_argument, 0,       'i'},
    {"mapping-method",       required_argument, 0,       'm'},
    {"mapping-options",      required_argument, 0,        2 },
    {"groups",               required_argument, 0,       'g'},
    {"task-allocations",     required_argument, 0,       't'},
    {"num-tasks",            required_argument, 0,        3 },
    {"num-task-allocations", required_argument, 0,        4 },
    {"check-accuracy",       no_argument,       0,        5 },
    {"realtime-clock",       no_argument,       0,        6 },
    {"verbose",              no_argument,       0,       'v'},
    {"show-gap-errors",      no_argument,       0,        7 },
    {nullptr,                0,                 nullptr,  0 }
  };

  ProfileOptions options;

  Stream groups_stream;
  Stream task_allocations_stream;

  for (;;) {
    int c = getopt_long(argc, argv, "hi:m:g:t:v", long_options, nullptr);
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
      case 'm':
        options.mapping_method.set(optarg);
        break;
      case 2:
        for (auto const &option : split(optarg, " "))
          options.mapping_options.set(option.c_str());
        break;
      case 'g':
        OPEN_STREAM(groups_stream, optarg);
        break;
      case 't':
        OPEN_STREAM(task_allocations_stream, optarg);
        break;
      case 3:
        options.num_tasks = stox<unsigned>(optarg);
        break;
      case 4:
        options.num_task_allocations = stox<unsigned>(optarg);
        break;
      case 5:
        options.check_accuracy = true;
        break;
      case 6:
        timer_realtime_enable();
        break;
      case 'v':
        ++options.verbosity;
        TIMER_ENABLE();
        break;
      case 7:
        options.show_gap_errors = true;
        break;
      default:
        return EXIT_FAILURE;
      }
    } catch (std::invalid_argument const &e) {
      error("invalid option argument:", e.what());
      return EXIT_FAILURE;
    }
  }

  CHECK_OPTION(options.library.is_set(), "--implementation option is mandatory");

  CHECK_OPTION(options.mapping_method.is_set(), "--mapping-method is mandatory");

  if (options.library.is("gap")) {
    CHECK_OPTION(!options.mapping_method.is("local_search"),
                 "local_search only supported when using mpsym");
  }

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

  if (options.check_accuracy)
    CHECK_OPTION(!options.library.is("gap"),
                 "--check-accuracy only available when using mpsym")

  try {
    profile(groups_stream, task_allocations_stream, options);
  } catch (std::exception const &e) {
    error("profiling failed:", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
