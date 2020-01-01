#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>

#include <getopt.h>
#include <libgen.h>

#include "arch_graph_system.h"
#include "arch_graph_cluster.h"

#include "profile_args.h"
#include "profile_parse.h"
#include "profile_read.h"
#include "profile_run.h"
#include "profile_util.h"

namespace
{

std::string progname;

void
usage(std::ostream &s)
{
  char const *opts[] = {
    "[-h|--help]",
    "-a|--arch-graphs ARCH_GRAPH",
    "[--disjoint-incomplete]",
    "[--disjoint-orbit-optimization]",
    "[--check-accuracy]",
    "[-v|--verbose]",
  };

  s << "usage: " << progname << '\n';
  for (char const *opt : opts)
    s << "  " << opt << '\n';
}

struct ProfileOptions
{
  bool disjoint_complete = true;
  bool disjoint_orbit_optimization = false;
  bool check_accuracy = false;
  bool verbose = false;
};

bool
is_cluster(std::shared_ptr<cgtl::ArchGraphSystem> ag)
{ return dynamic_cast<cgtl::ArchGraphCluster *>(ag.get()); }

std::shared_ptr<cgtl::ArchGraphCluster>
as_cluster(std::shared_ptr<cgtl::ArchGraphSystem> ag)
{ return std::dynamic_pointer_cast<cgtl::ArchGraphCluster>(ag); }

std::vector<cgtl::PermGroup>
decompose_cluster(std::shared_ptr<cgtl::ArchGraphCluster> ag,
                  ProfileOptions const &options)
{
  return ag->automorphisms().disjoint_decomposition(
    options.disjoint_complete, options.disjoint_orbit_optimization);
}

void
decompose_cluster_wrapper(std::shared_ptr<cgtl::ArchGraphCluster> ag,
                          ProfileOptions const &options,
                          double *t)
{
  using cgtl::PermGroup;

  auto decomposition(
    run_cpp([&]{ return decompose_cluster(ag, options); }, t));

  if (options.verbose)
    info("Decomposes into", decomposition.size(), "clusters");

  if (decomposition.size() < ag->num_subsystems())
    warning("Expected decomposition into at least", ag->num_subsystems(),
            "subsystems but found only", decomposition.size());

  if (options.check_accuracy) {
    info("Checking accuracy...");

    auto reconstruction(
      PermGroup::group_union(decomposition.begin(), decomposition.end()));

    if (reconstruction != ag->automorphisms())
      info("Decomposition is incorrect");
    else
      info("Decomposition is correct");
  }
}

double
run(std::shared_ptr<cgtl::ArchGraphSystem> ag,
    ProfileOptions const &options)
{
  double t;

  if (is_cluster(ag))
    decompose_cluster_wrapper(as_cluster(ag), options, &t);

  return t;
}

void
profile(Stream &arch_graphs_stream,
        ProfileOptions const &options)
{
  foreach_line(arch_graphs_stream.stream,
               [&](std::string const &line, unsigned lineno){

    auto ag(parse_arch_graph_system(line));

    if (options.verbose) {
      info("Decomposing arch graph", lineno);
      info("==>", ag->num_processors(), "processors");
      if (is_cluster(ag))
        info("==>", as_cluster(ag)->num_subsystems(), "subsystems");
    } else {
      info("Decomposing arch graph", lineno);
    }

    double t = run(ag, options);

    result("Runtime:", t, "s");
  });
}

} // namespace


int
main(int argc, char **argv)
{
  progname = basename(argv[0]);

  struct option long_options[] = {
    {"help",                        no_argument,       0,       'h'},
    {"arch_graphs",                 no_argument,       0,       'a'},
    {"disjoint-incomplete",         no_argument,       0,        1 },
    {"disjoint-orbit-optimization", no_argument,       0,        2 },
    {"check-accuracy",              no_argument,       0,        3 },
    {"verbose",                     no_argument,       0,       'v'},
    {nullptr,                       0,                 nullptr,  0 }
  };

  ProfileOptions options;

  Stream arch_graphs_stream;

  for (;;) {
    int c = getopt_long(argc, argv, "ha:v", long_options, nullptr);
    if (c == -1)
      break;

    try {
      switch (c) {
      case 'h':
        usage(std::cout);
        return EXIT_SUCCESS;
      case 'a':
        OPEN_STREAM(arch_graphs_stream, optarg);
        break;
      case 1:
        options.disjoint_complete = false;
        break;
      case 2:
        options.disjoint_orbit_optimization = true;
        break;
      case 3:
        options.check_accuracy = true;
        break;
      case 'v':
        options.verbose = true;
        break;
      default:
        return EXIT_FAILURE;
      }
    } catch (std::invalid_argument const &e) {
      error("invalid option argument:", e.what());
      return EXIT_FAILURE;
    }
  }

  profile(arch_graphs_stream, options);

  return EXIT_SUCCESS;
}
