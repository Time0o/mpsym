#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <functional>
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
#include "task_mapping.h"
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
    "-m|--repr-method {iterate|orbits|local_search}",
    "--repr-variant {local_search_bfs|local_search_dfs|local_search_sa_linear}",
    "--repr-local-search-invert-generators",
    "--repr-local-search-append-generators",
    "--repr-local-search-iterations",
    "--repr-local-search-sa-T-init",
    "[--repr-options {dont_decompose,dont_match,dont_optimize_symmetric}]",
    "[-g|--groups GROUPS]",
    "[-a|--arch-graph ARCH_GRAPH]",
    "[--arch-graph-args ARCH_GRAPH_ARGS]",
    "[--dont-decompose-arch-graph]",
    "-t|--task-mappings TASK_ALLOCATIONS",
    "[-l|--task-mappings-limit TASK_ALLOCATIONS_LIMIT]",
    "[-r|--num-runs NUM_RUNS]",
    "[--num-discarded-runs NUM_DISCARDED_RUNS]",
    "[--summarize-runs]",
    "[-c|--check-accuracy-gap]",
    "[--check-accuracy-mpsym]",
    "[-v|--verbose]",
    "[--compile-gap]",
    "[--show-gap-errors]"
  };

  s << "usage: " << progname << '\n';
  for (char const *opt : opts)
    s << "  " << opt << '\n';
}

struct ProfileOptions
{
  VariantOption library{"gap", "mpsym"};
  VariantOption repr_method{"iterate", "orbits", "local_search"};
  VariantOption repr_variant{
    "local_search_bfs", "local_search_dfs", "local_search_sa_linear"};
  VariantOptionSet repr_options{
    "dont_decompose", "dont_match", "dont_optimize_symmetric"};

  bool repr_local_search_invert_generators = false;
  unsigned repr_local_search_append_generators = 0u;
  double repr_local_search_sa_iterations = 0.0;
  double repr_local_search_sa_T_init = 0.0;

  bool groups_input = false;
  bool arch_graph_input = false;
  std::vector<std::string> arch_graph_args;
  bool dont_decompose_arch_graph = false;
  unsigned task_mappings_limit = 0u;
  unsigned num_runs = 1u;
  unsigned num_discarded_runs = 0u;
  bool summarize_runs = false;
  bool check_accuracy_gap = false;
  bool check_accuracy_mpsym = false;
  int verbosity = 0;
  bool compile_gap = false;
  bool show_gap_errors = false;
};

std::string map_tasks_gap_local_search(ProfileOptions const &options)
{
  return R"(
orbit_repr:=task_mapping;

generators:=GeneratorsOfGroup(automorphisms);
possible_representatives:=EmptyPlist(Size(generators));

for i in [1..Size(generators)] do
  stationary:=true;

  permuted:=OnTuples(orbit_repr, generators[i]);

  if permuted < orbit_repr then
    possible_representatives[i]:=OnTuples(orbit_repr, generators[i]);
    stationary:=false;
  fi;

  if stationary then
    break;
  fi;

  orbit_repr:=Minimum(possible_representatives);
od;

if HTAdd(orbit_representatives_hash, orbit_repr, true) <> fail then
  Append(orbit_representatives, [orbit_repr]);
fi;
  )";
}

std::string map_tasks_gap_iterate(ProfileOptions const &options)
{
  std::string ret(R"(
orbit_repr:=task_mapping;
orbit_repr_new:=true;

for element in automorphisms do
  permuted:=OnTuples(task_mapping, element);
  )");

  if (options.repr_options.is_set("dont_match")) {
    ret += R"(
  if permuted < orbit_repr then
    orbit_repr:=permuted;
  fi;
od;

if HTAdd(orbit_representatives_hash, orbit_repr, true) <> fail then
  Append(orbit_representatives, [orbit_repr]);
fi;
    )";
  } else {
    ret += R"(
  if HTValue(orbit_representatives_hash, permuted) <> fail then
    orbit_repr_new:=false;
    break;
  elif permuted < orbit_repr then
    orbit_repr:=permuted;
  fi;
od;

if orbit_repr_new then
  HTAdd(orbit_representatives_hash, orbit_repr, true);
  Append(orbit_representatives, [orbit_repr]);
fi;
    )";
  }

  return ret;
}

std::string map_tasks_gap_orbits(ProfileOptions const &options)
{
  std::string ret;

  if (options.repr_options.is_set("dont_match")) {
    ret = R"(
orbit:=Orb(automorphisms, task_mapping, OnTuples);
orbit_repr:=Elements(Enumerate(orbit))[1];
    )";
  } else {
    ret = R"(
orbit_options:=rec(lookingfor:=orbit_representatives_hash);
orbit:=Orb(automorphisms, task_mapping, OnTuples, orbit_options);
orbit_repr:=Elements(Enumerate(orbit))[1];
    )";
  }

  ret += R"(
if HTAdd(orbit_representatives_hash, orbit_repr, true) <> fail then
  Append(orbit_representatives, [orbit_repr]);
fi;
  )";

  return ret;
}

std::string map_tasks_gap(ProfileOptions const &options)
{
  std::stringstream ss;

  if (options.verbosity > 0)
    ss << "Print(\"DEBUG: Constructing BSGS\\n\");\n";

  ss << "StabChain(automorphisms);\n";

  ss << "orbit_representatives:=[];\n";
  ss << "orbit_representatives_hash:=HTCreate([1,2,3]);\n";

  // map tasks mappings one by one
  ss << "n:=1;\n";
  ss << "for task_mapping in task_mappings do\n";

  // display progress
  if (options.verbosity > 0) {
    ss << "    Print(\"DEBUG: Mapping task \", n, \" of \", "
       << "          Length(task_mappings), \"\\r\\c\");\n";
  }

  // concrete code depending on chosen implementation
  if (options.repr_method.is("iterate")) {
    ss << map_tasks_gap_iterate(options);
  } else if (options.repr_method.is("orbits")) {
    ss << map_tasks_gap_orbits(options);
  } else if (options.repr_method.is("local_search")) {
    ss << map_tasks_gap_local_search(options);
  } else {
    throw std::logic_error("unreachable");
  }

  ss << "  n:=n+1;\n";
  ss << "od;\n";

  // display orbit representatives found
  if (options.check_accuracy_gap || options.verbosity > 0) {
    ss << "Print(\"\\n;DEBUG: => Found \", Length(orbit_representatives), "
       << "      \" orbit representatives;\\n\");\n";

    if (options.check_accuracy_gap || options.verbosity > 1) {
      ss << "for orbit_repr in orbit_representatives do\n";
      ss << "  Print(\"DEBUG: \", orbit_repr, \";\\n\");\n";
      ss << "od;\n";
    }
  }

  return ss.str();
}

mpsym::ArchGraphSystem::ReprOptions map_tasks_mpsym_repr_options(
  ProfileOptions const &options)
{
  using mpsym::ArchGraphSystem;

  ArchGraphSystem::ReprOptions repr_options;

  if (options.repr_method.is("iterate")) {
    repr_options.method = ArchGraphSystem::ReprMethod::ITERATE;
  } else if (options.repr_method.is("orbits")) {
    repr_options.method = ArchGraphSystem::ReprMethod::ORBITS;
  } else if (options.repr_method.is("local_search")) {
    repr_options.method = ArchGraphSystem::ReprMethod::LOCAL_SEARCH;

    if (options.repr_variant.is("local_search_bfs"))
      repr_options.variant = ArchGraphSystem::ReprVariant::LOCAL_SEARCH_BFS;
    else if (options.repr_variant.is("local_search_dfs"))
      repr_options.variant = ArchGraphSystem::ReprVariant::LOCAL_SEARCH_DFS;
    else if (options.repr_variant.is("local_search_sa_linear"))
      repr_options.variant = ArchGraphSystem::ReprVariant::LOCAL_SEARCH_SA_LINEAR;

    repr_options.local_search_invert_generators =
      options.repr_local_search_invert_generators;

    repr_options.local_search_append_generators =
      options.repr_local_search_append_generators;

    if (options.repr_local_search_sa_iterations > 0u) {
      repr_options.local_search_sa_iterations =
        options.repr_local_search_sa_iterations;
    }

    if (options.repr_local_search_sa_T_init > 0.0) {
      repr_options.local_search_sa_T_init =
        options.repr_local_search_sa_T_init;
    }

  } else {
    throw std::logic_error("unreachable");
  }

  if (options.repr_options.is_set("dont_match"))
    repr_options.match = false;

  if (options.repr_options.is_set("dont_optimize_symmetric"))
    repr_options.optimize_symmetric = false;

  return repr_options;
}

mpsym::TaskOrbits map_tasks_mpsym(
  std::shared_ptr<mpsym::ArchGraphSystem> ags,
  mpsym::TaskMappingVector const &task_mappings,
  mpsym::ArchGraphSystem::ReprOptions const &repr_options,
  ProfileOptions const &options)
{
  using mpsym::ArchGraphAutomorphisms;
  using mpsym::ArchGraphSystem;
  using mpsym::PermGroup;
  using mpsym::TaskOrbits;

  if (options.verbosity > 0)
    debug("Constructing BSGS");

  ags->init_repr();

  TaskOrbits task_orbits;
  for (auto i = 0u; i < task_mappings.size(); ++i) {
    if (options.verbosity > 0)
      debug_progress("Mapping task", i + 1u, "of", task_mappings.size());

    ags->repr(task_mappings[i], &task_orbits, &repr_options);
  }

  return task_orbits;
}

void map_tasks_gap_wrapper(gap::PermGroup const &automorphisms,
                           gap::TaskMappingVector const &task_mappings,
                           ProfileOptions const &options,
                           mpsym::TaskOrbits *task_orbits,
                           std::vector<double> *ts)
{
  using mpsym::TaskOrbits;

  // run gap script

  auto gap_automorphisms("automorphisms:=" + automorphisms + ";");
  auto gap_task_mappings("task_mappings:=[\n" + task_mappings + "\n];");
  auto gap_script(map_tasks_gap(options));

  auto gap_output(run_gap({"orb", "grape"},
                          {{"automorphisms", gap_automorphisms, options.compile_gap},
                           {"task_mappings", gap_task_mappings, false}},
                          gap_script,
                          options.num_discarded_runs,
                          options.num_runs,
                          options.verbosity == 0,
                          !options.show_gap_errors,
                          options.compile_gap,
                          ts));

  // parse output

  if (options.check_accuracy_gap) {
    auto representatives_gap(parse_task_mappings_gap_to_mpsym(
      std::vector<std::string>(gap_output.begin() + 2, gap_output.end())));

    task_orbits->insert_all(representatives_gap.begin(),
                            representatives_gap.end());
  }
}

void map_tasks_mpsym_wrapper(std::shared_ptr<mpsym::ArchGraphSystem> ags,
                             mpsym::TaskMappingVector const &task_mappings,
                             ProfileOptions const &options,
                             mpsym::TaskOrbits *task_orbits,
                             std::vector<double> *ts)
{
  using mpsym::TaskOrbits;

  auto repr_options(map_tasks_mpsym_repr_options(options));

  *task_orbits = run_cpp(
    [&]{ return map_tasks_mpsym(ags, task_mappings, repr_options, options); },
    options.num_discarded_runs,
    options.num_runs,
    ts);

  if (options.verbosity > 0) {
    debug_progress_done();

    debug("=> Found", task_orbits->num_orbits(), "orbit representatives");
    if (options.verbosity > 1) {
      for (auto const &repr : *task_orbits)
        debug(DUMP(repr));
    }

    std::cout << std::endl;
  }
}

void check_accuracy(mpsym::TaskOrbits const &task_orbits_actual,
                    mpsym::TaskOrbits const &task_orbits_check,
                    ProfileOptions const &options)
{
  using mpsym::TaskMapping;

  if (task_orbits_actual == task_orbits_check) {
    info("Orbit representatives match");
    return;
  }

  info("Orbit representatives do not match:");

  // construct representative sets
  std::set<TaskMapping> reprs_actual, reprs_check,
                        reprs_found, reprs_missing, reprs_extra;

  for (auto const &repr : task_orbits_actual)
    reprs_actual.insert(repr);

  for (auto const &repr : task_orbits_check)
    reprs_check.insert(repr);

  // find missing/extra representatives
  if (!reprs_check.empty()) {
    for (auto const &repr : reprs_check) {
      if (reprs_actual.find(repr) == reprs_actual.end())
        reprs_missing.insert(repr);
      else
        reprs_found.insert(repr);
    }

    info("=>", reprs_missing.size(), "Missing orbit representatives");
    if (options.verbosity > 1) {
      for (auto const &repr : reprs_missing)
        info(DUMP(repr));
    }
  }

  if (!reprs_actual.empty()) {
    for (auto const &repr : reprs_actual) {
      if (reprs_check.find(repr) == reprs_check.end())
        reprs_extra.insert(repr);
    }

    info("=>", reprs_extra.size(), "Extra orbit representatives");
    if (options.verbosity > 1) {
      for (auto const &repr : reprs_extra)
        info(DUMP(repr));
    }
  }

  info("=> Found", reprs_found.size(), "of", reprs_check.size(), "representatives");
}

void run(std::shared_ptr<mpsym::ArchGraphSystem> ags,
         std::shared_ptr<mpsym::ArchGraphSystem> ags_check,
         std::string const &task_mappings,
         ProfileOptions const &options)
{
  using namespace std::placeholders;

  using mpsym::TaskOrbits;

  std::vector<double> ts;

  auto run_gap(std::bind(map_tasks_gap_wrapper,
                         _1,
                         parse_task_mappings_gap(task_mappings),
                         _2,
                         _3,
                         _4));

  auto run_mpsym(std::bind(map_tasks_mpsym_wrapper,
                           _1,
                           parse_task_mappings_mpsym(task_mappings),
                           _2,
                           _3,
                           _4));

  if (options.library.is("gap")) {
    run_gap(ags->to_gap(), options, nullptr, &ts);

  } else if (options.library.is("mpsym")) {
    TaskOrbits task_orbits, task_orbits_check;

    run_mpsym(ags, options, &task_orbits, &ts);

    if (options.check_accuracy_gap || options.check_accuracy_mpsym) {
      info("Checking accuracy...");

      auto options_(options);
      options_.repr_method.set("orbits");
      options_.repr_options.unset("dont_match");
      options_.repr_options.unset("dont_optimize_symmetric");
      options_.verbosity = 0;

      if (options.check_accuracy_gap)
        run_gap(ags_check->to_gap(), options_, &task_orbits_check, nullptr);
      else if (options.check_accuracy_mpsym)
        run_mpsym(ags_check, options_, &task_orbits_check, nullptr);

      check_accuracy(task_orbits, task_orbits_check, options);
    }
  }

  dump_runs(ts, options.summarize_runs);
}

void do_profile(Stream &automorphisms_stream,
                Stream &task_mappings_stream,
                ProfileOptions const &options)
{
  using mpsym::ArchGraphAutomorphisms;
  using mpsym::ArchGraphSystem;
  using mpsym::PermGroup;

  std::shared_ptr<ArchGraphSystem> ags, ags_check;

  auto task_mappings(read_file(task_mappings_stream.stream,
                               options.task_mappings_limit));

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
      ags_check = ags;
    });

  } else if (options.arch_graph_input) {
    auto arch_graph(read_file(automorphisms_stream.stream));

    ags = ArchGraphSystem::from_lua(arch_graph, options.arch_graph_args);
    ags_check = ags;

    if (options.dont_decompose_arch_graph ||
        options.repr_options.is_set("dont_decompose")) {

      if (options.verbosity > 0)
        debug("Determining automorphisms");

      ags = std::make_shared<ArchGraphAutomorphisms>(ags->automorphisms());
    }
  }

  run(ags, ags_check, task_mappings, options);
}

} // namespace

int main(int argc, char **argv)
{
  progname = basename(argv[0]);

  struct option long_options[] = {
    {"help",                                no_argument,       0,       'h'},
    {"implementation",                      required_argument, 0,       'i'},
    {"repr-method",                         required_argument, 0,       'm'},
    {"repr-variant",                        required_argument, 0,        1 },
    {"repr-local-search-invert-generators", no_argument,       0,        2 },
    {"repr-local-search-append-generators", required_argument, 0,        3 },
    {"repr-local-search-sa-iterations",     required_argument, 0,        4 },
    {"repr-local-search-sa-T-init",         required_argument, 0,        5 },
    {"repr-options",                        required_argument, 0,        6 },
    {"groups",                              required_argument, 0,       'g'},
    {"arch-graph",                          required_argument, 0,       'a'},
    {"arch-graph-args",                     required_argument, 0,        7 },
    {"dont-decompose-arch-graph",           no_argument, 0,              8 },
    {"task-mappings",                       required_argument, 0,       't'},
    {"task-mappings-limit",                 required_argument, 0,       'l'},
    {"num-runs",                            required_argument, 0,       'r'},
    {"num-discarded-runs",                  required_argument, 0,        9 },
    {"summarize-runs",                      required_argument, 0,        10},
    {"check-accuracy-gap",                  no_argument,       0,       'c'},
    {"check-accuracy-mpsym",                no_argument,       0,        11},
    {"verbose",                             no_argument,       0,       'v'},
    {"compile-gap",                         no_argument,       0,        12},
    {"show-gap-errors",                     no_argument,       0,        13},
    {nullptr,                               0,                 nullptr,  0 }
  };

  ProfileOptions options;

  Stream automorphisms_stream;
  Stream task_mappings_stream;

  for (;;) {
    int c = getopt_long(argc, argv, "hi:m:g:a:t:l:r:cv", long_options, nullptr);
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
      case 1:
        options.repr_variant.set(optarg);
        break;
      case 2:
        options.repr_local_search_invert_generators = true;
        break;
      case 3:
        options.repr_local_search_append_generators = stox<unsigned>(optarg);
        break;
      case 4:
        options.repr_local_search_sa_iterations = stox<unsigned>(optarg);
        break;
      case 5:
        options.repr_local_search_sa_T_init = stof<double>(optarg);
        break;
      case 6:
        foreach_option(optarg,
                       [&](std::string const &option)
                       { options.repr_options.set(option.c_str()); });
        break;
      case 'g':
        OPEN_STREAM(automorphisms_stream, optarg);
        options.groups_input = true;
        break;
      case 'a':
        OPEN_STREAM(automorphisms_stream, optarg);
        options.arch_graph_input = true;
        break;
      case 7:
        foreach_option(optarg,
                       [&](std::string const &option)
                       { options.arch_graph_args.push_back(option); });
        break;
      case 8:
        options.dont_decompose_arch_graph = true;
        break;
      case 't':
        OPEN_STREAM(task_mappings_stream, optarg);
        break;
      case 'l':
        options.task_mappings_limit = stox<unsigned>(optarg);
        break;
      case 'r':
        options.num_runs = stox<unsigned>(optarg);
        break;
      case 9:
        options.num_discarded_runs = stox<unsigned>(optarg);
        break;
      case 10:
        options.summarize_runs = true;
        break;
      case 'c':
        options.check_accuracy_gap = true;
        break;
      case 11:
        options.check_accuracy_mpsym = true;
        break;
      case 'v':
        ++options.verbosity;
        TIMER_ENABLE();
        break;
      case 12:
        options.compile_gap = true;
        break;
      case 13:
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

  CHECK_OPTION(task_mappings_stream.valid,
               "--task-mappings option is mandatory");

  CHECK_OPTION(options.groups_input != options.arch_graph_input,
               "EITHER --arch-graph OR --groups must be given");

  CHECK_OPTION(!options.library.is("gap") ||
               !(options.check_accuracy_gap || options.check_accuracy_mpsym),
               "--check-accuracy-* only available when using mpsym");

  try {
    do_profile(automorphisms_stream, task_mappings_stream, options);
  } catch (std::exception const &e) {
    error("profiling failed:", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
