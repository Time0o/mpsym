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

#include "arch_graph_automorphisms.h"
#include "arch_graph_system.h"
#include "dump.h"
#include "perm_set.h"
#include "task_allocation.h"
#include "task_orbits.h"
#include "timer.h"

#include "profile_args.h"
#include "profile_parse.h"
#include "profile_read.h"
#include "profile_run.h"
#include "profile_timer.h"
#include "profile_util.h"

using namespace profile;

namespace
{

std::string progname;

void usage(std::ostream &s)
{
  char const *opts[] = {
    "[-h|--help]",
    "-i|--implementation {gap|mpsym}",
    "-m|--repr-method {iterate|local_search|orbits}",
    "[--repr-options {dont_match_reprs}]",
    "[-g|--groups GROUPS]",
    "[-a|--arch-graph ARCH_GRAPH]",
    "-t|--task-allocations TASK_ALLOCATIONS",
    "[-l|--task-allocations-limit TASK_ALLOCATIONS_LIMIT]",
    "[-c|--check-accuracy]",
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
  VariantOption repr_method{"iterate", "local_search", "orbits"};
  VariantOptionSet repr_options{"dont_match_reprs"};
  bool groups_input = false;
  bool arch_graph_input = false;
  unsigned task_allocation_limit = 0u;
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

  if (options.repr_options.is_set("dont_match_reprs")) {
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

  if (options.repr_options.is_set("dont_match_reprs")) {
    ss << "orbit:=Orb(automorphisms, task_allocation, OnTuples);\n";
    ss << "orbit_repr:=Elements(Enumerate(orbit))[1];\n";

  } else {
    ss << "orbit_options:=rec(lookingfor:=orbit_representatives_hash);\n";
    ss << "orbit:=Orb(automorphisms, task_allocation, OnTuples, orbit_options);\n";
    ss << "orbit_repr:=Elements(Enumerate(orbit))[1];\n";
  }

  ss << "if HTAdd(orbit_representatives_hash, orbit_repr, true) <> fail then\n";
  ss << "  Append(orbit_representatives, [orbit_repr]);\n";
  ss << "fi;\n";

  return ss.str();
}

std::string map_tasks_gap(
  gap::PermGroup const &automorphisms,
  gap::TaskAllocationVector const &task_allocations,
  ProfileOptions const &options)
{
  std::stringstream ss;

  // construct the automorphism group
  ss << "automorphisms:=" << automorphisms << ";\n";

  // construct the vector of task allocations to be mapped
  ss << "task_allocations:=[\n";
  ss << task_allocations;
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

  // concrete code depending on chosen implementation
  if (options.repr_method.is("iterate"))
    ss << map_tasks_gap_iterate(options);
  else if (options.repr_method.is("orbits"))
    ss << map_tasks_gap_orbits(options);
  else
    throw std::logic_error("unreachable");

  ss << "  n:=n+1;\n";
  ss << "od;\n";

  // display orbit representatives found
  if (options.check_accuracy || options.verbosity > 0) {
    ss << "Print(\"\\n;DEBUG: => Found \", Length(orbit_representatives), "
          "      \" orbit representatives;\\n\");\n";

    if (options.check_accuracy || options.verbosity > 1) {
      ss << "for orbit_repr in orbit_representatives do\n";
      ss << "  Print(\"DEBUG: \", orbit_repr, \";\\n\");\n";
      ss << "od;\n";
    }
  }

  return ss.str();
}

mpsym::TaskOrbits map_tasks_mpsym(
  std::shared_ptr<mpsym::ArchGraphSystem> ags,
  mpsym::TaskAllocationVector const &task_allocations,
  ProfileOptions const &options)
{
  using mpsym::ArchGraphAutomorphisms;
  using mpsym::ArchGraphSystem;
  using mpsym::PermGroup;
  using mpsym::TaskOrbits;

  ArchGraphSystem::ReprOptions repr_options;

  if (options.repr_method.is("iterate"))
    repr_options.method = ArchGraphSystem::ReprMethod::ITERATE;
  else if (options.repr_method.is("local_search"))
    repr_options.method = ArchGraphSystem::ReprMethod::LOCAL_SEARCH;
  else if (options.repr_method.is("orbits"))
    repr_options.method = ArchGraphSystem::ReprMethod::ORBITS;
  else
    throw std::logic_error("unreachable");

  if (options.repr_options.is_set("dont_match_reprs"))
      repr_options.match_reprs = false;

  TaskOrbits task_orbits;
  for (auto i = 0u; i < task_allocations.size(); ++i) {
    if (options.verbosity > 0)
      debug_progress("Mapping task", i + 1u, "of", task_allocations.size());

    ags->repr(task_allocations[i], &repr_options, &task_orbits);
  }

  return task_orbits;
}

void map_tasks_gap_wrapper(gap::PermGroup const &automorphisms,
                           gap::TaskAllocationVector const &task_allocations,
                           ProfileOptions const &options,
                           mpsym::TaskOrbits *task_orbits,
                           double *t)
{
  using mpsym::TaskOrbits;

  // run gap script

  auto gap_script(map_tasks_gap(automorphisms,
                                task_allocations,
                                options));

  std::vector<double> ts;

  auto gap_output(run_gap({"orb", "grape"},
                          gap_script,
                          0,
                          1,
                          options.check_accuracy,
                          !options.show_gap_errors,
                          &ts));

  if (t)
    *t = ts[0];

  // parse output

  if (options.check_accuracy) {
    auto representatives_gap(parse_task_allocations_gap_to_mpsym(
      std::vector<std::string>(gap_output.begin() + 2, gap_output.end())));

    task_orbits->insert_all(representatives_gap.begin(),
                            representatives_gap.end());
  }
}

void map_tasks_mpsym_wrapper(std::shared_ptr<mpsym::ArchGraphSystem> ags,
                             mpsym::TaskAllocationVector const &task_allocations,
                             ProfileOptions const &options,
                             mpsym::TaskOrbits *task_orbits,
                             double *t)
{
  using mpsym::TaskOrbits;

  std::vector<double> ts;

  *task_orbits = run_cpp([&]{
    return map_tasks_mpsym(ags, task_allocations, options);
  }, 0, 1, &ts);

  if (t)
    *t = ts[0];

  if (options.verbosity > 0) {
    debug_progress_done();

    debug("=> Found", task_orbits->num_orbits(), "orbit representatives");
    if (options.verbosity > 1) {
      for (auto const &repr : *task_orbits)
        debug(DUMP(repr));
    }
  }
}

void check_accuracy(mpsym::TaskOrbits const &task_orbits_mpsym,
                    mpsym::TaskOrbits const &task_orbits_gap,
                    ProfileOptions const &options)
{
  using mpsym::TaskAllocation;

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

double run(std::shared_ptr<mpsym::ArchGraphSystem> ags,
           std::string const &task_allocations,
           ProfileOptions const &options)
{
  using mpsym::TaskOrbits;

  double t;

  if (options.library.is("gap")) {
    map_tasks_gap_wrapper(ags->to_gap(),
                          parse_task_allocations_gap(task_allocations),
                          options,
                          nullptr,
                          &t);

  } else if (options.library.is("mpsym")) {
    TaskOrbits task_orbits_mpsym, task_orbits_gap;

    map_tasks_mpsym_wrapper(ags,
                            parse_task_allocations_mpsym(task_allocations),
                            options,
                            &task_orbits_mpsym,
                            &t);

    if (options.check_accuracy) {
      info("Checking accuracy...");

      map_tasks_gap_wrapper(ags->to_gap(),
                            parse_task_allocations_gap(task_allocations),
                            options,
                            &task_orbits_gap,
                            nullptr);

      check_accuracy(task_orbits_mpsym, task_orbits_gap, options);
    }
  }

  return t;
}

void do_profile(Stream &automorphisms_stream,
                Stream &task_allocations_stream,
                ProfileOptions const &options)
{
  using mpsym::ArchGraphAutomorphisms;
  using mpsym::ArchGraphSystem;
  using mpsym::PermGroup;

  std::shared_ptr<ArchGraphSystem> ags;

  auto task_allocations(read_file(task_allocations_stream.stream,
                                  options.task_allocation_limit));

  if (options.verbosity > 0)
    debug("Implementation:", options.library.get());

  if (options.groups_input) {
    foreach_line(automorphisms_stream.stream,
                 [&](std::string const &line, unsigned lineno){

      auto group(parse_group(line));

      if (options.verbosity > 0) {
        info("Using automorphism group", lineno);
        info("=> degree", group.degree);
        info("=> orders", group.order);
        info("=> generators", group.generators);
      } else {
        info("Using automorphism group", lineno);
      }

      ags = group.to_arch_graph_system();
    });

  } else if (options.arch_graph_input) {
    ags = ArchGraphSystem::from_lua(read_file(automorphisms_stream.stream));
  }

  double t = run(ags, task_allocations, options);

  result("Runtime:", t, "s");

  if (options.verbosity > 0) {
    debug("Timer dumps:");
    if (options.repr_method.is("iterate"))
      debug_timer_dump("map bruteforce iterate");
    else if (options.repr_method.is("local_search"))
      debug_timer_dump("map approx local search");
    else if (options.repr_method.is("orbits"))
      debug_timer_dump("map bruteforce orbits");
  }
}

} // namespace

int main(int argc, char **argv)
{
  progname = basename(argv[0]);

  struct option long_options[] = {
    {"help",                   no_argument,       0,       'h'},
    {"implementation",         required_argument, 0,       'i'},
    {"repr-method",            required_argument, 0,       'm'},
    {"repr-options",           required_argument, 0,        2 },
    {"groups",                 required_argument, 0,       'g'},
    {"arch-graph",             required_argument, 0,       'a'},
    {"task-allocations",       required_argument, 0,       't'},
    {"task-allocations-limit", required_argument, 0,       'l'},
    {"check-accuracy",         no_argument,       0,       'c'},
    {"verbose",                no_argument,       0,       'v'},
    {"show-gap-errors",        no_argument,       0,        3 },
    {nullptr,                  0,                 nullptr,  0 }
  };

  ProfileOptions options;

  Stream automorphisms_stream;
  Stream task_allocations_stream;

  for (;;) {
    int c = getopt_long(argc, argv, "hi:m:g:a:t:l:v", long_options, nullptr);
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
        options.repr_method.set(optarg);
        break;
      case 2:
        for (auto const &option : split(optarg, " "))
          options.repr_options.set(option.c_str());
        break;
      case 'g':
        OPEN_STREAM(automorphisms_stream, optarg);
        options.groups_input = true;
        break;
      case 'a':
        OPEN_STREAM(automorphisms_stream, optarg);
        options.arch_graph_input = true;
        break;
      case 't':
        OPEN_STREAM(task_allocations_stream, optarg);
        break;
      case 'l':
        options.task_allocation_limit = stox<unsigned>(optarg);
        break;
      case 'c':
        options.check_accuracy = true;
        break;
      case 'v':
        ++options.verbosity;
        TIMER_ENABLE();
        break;
      case 3:
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

  CHECK_OPTION(options.repr_method.is_set(), "--repr-method is mandatory");

  CHECK_OPTION(!options.library.is("gap") || !options.repr_method.is("local_search"),
               "local_search only supported when using mpsym");

  CHECK_OPTION(task_allocations_stream.valid,
               "--task-allocations option is mandatory");

  CHECK_OPTION(options.groups_input != options.arch_graph_input,
               "EITHER --arch-graph OR --groups must be given");

  CHECK_OPTION(!options.check_accuracy || !options.library.is("gap"),
               "--check-accuracy only available when using mpsym")

  try {
    do_profile(automorphisms_stream, task_allocations_stream, options);
  } catch (std::exception const &e) {
    error("profiling failed:", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
