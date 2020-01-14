#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <unordered_set>

#include <getopt.h>
#include <libgen.h>

#include "arch_graph_system.h"
#include "arch_graph_cluster.h"
#include "arch_uniform_super_graph.h"
#include "cartesian_product.h"

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

bool
is_supergraph(std::shared_ptr<cgtl::ArchGraphSystem> ag)
{ return dynamic_cast<cgtl::ArchUniformSuperGraph *>(ag.get()); }

std::shared_ptr<cgtl::ArchUniformSuperGraph>
as_supergraph(std::shared_ptr<cgtl::ArchGraphSystem> ag)
{ return std::dynamic_pointer_cast<cgtl::ArchUniformSuperGraph>(ag); }

std::vector<cgtl::PermGroup>
decompose_cluster(std::shared_ptr<cgtl::ArchGraphCluster> ag,
                  ProfileOptions const &options)
{
  return ag->automorphisms().disjoint_decomposition(
    options.disjoint_complete, options.disjoint_orbit_optimization);
}

std::vector<cgtl::PermGroup>
decompose_supergraph(std::shared_ptr<cgtl::ArchUniformSuperGraph> ag,
                     ProfileOptions const &)
{ return ag->automorphisms().wreath_decomposition(); }

void
decompose_cluster_wrapper(std::shared_ptr<cgtl::ArchGraphCluster> ag,
                          ProfileOptions const &options,
                          double *t)
{
  using cgtl::PermGroup;

  if (options.verbose)
    info("Trying to decompose cluster...");

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

void
decompose_supergraph_wrapper(std::shared_ptr<cgtl::ArchUniformSuperGraph> ag,
                             ProfileOptions const &options,
                             double *t)
{
  using cgtl::Perm;
  using cgtl::PermGroup;

  if (options.verbose)
    info("Trying to decompose supergraph...");

  auto decomposition(
    run_cpp([&]{ return decompose_supergraph(ag, options); }, t));

  if (decomposition.empty())
    warning("Failed to find supergraph decomposition");

  if (options.check_accuracy) {
    info("Checking accuracy...");

    // determine all elements in group reconstructed from wreath decomposition
    std::unordered_set<Perm> decomposition_elements;

    std::vector<std::vector<Perm>> sigmas;

    for (PermGroup const &pg : decomposition) {
      std::vector<Perm> sigma;

      // TODO
      for (Perm const &perm : pg)
        sigma.push_back(perm);

      sigmas.push_back(sigma);
    }

    std::vector<Perm> chain;
    while (next_in_cartesian_product(sigmas.begin(), sigmas.end(), chain.begin())) {
      Perm perm(chain[0]);
      for (auto i = 1u; i < chain.size(); ++i)
        perm *= chain[1];

      decomposition_elements.insert(perm);
    }

    // order sanity check
    if (decomposition_elements.size() != ag->automorphisms().order()) {
      // determine all elements in automorphism group
      std::unordered_set<Perm> automorphisms_elements;

      for (Perm const &perm : ag->automorphisms())
        automorphisms_elements.insert(perm);

      if (automorphisms_elements != decomposition_elements)
        info("Decomposition is incorrect, elements do not match");

    } else {
      info("Decomposition is incorrect, expected", ag->automorphisms().order(),
            "elements but got", decomposition_elements.size());
    }
  }
}

double
run(std::shared_ptr<cgtl::ArchGraphSystem> ag,
    ProfileOptions const &options)
{
  double t;

  if (is_cluster(ag)) {
    if (options.verbose)
      info("Graph is cluster");

    decompose_cluster_wrapper(as_cluster(ag), options, &t);

  } else if (is_supergraph(ag)) {
    if (options.verbose)
      info("Graph is supergraph");

    decompose_supergraph_wrapper(as_supergraph(ag), options, &t);
  }

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
      info("Decomposing graph", lineno);
      info("=>", ag->num_processors(), "processors");
    } else {
      info("Decomposing graph", lineno);
    }

    if (options.verbose)
      info("Constructing automorphism group");

    (void)ag->automorphisms();

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

  CHECK_OPTION(arch_graphs_stream.valid, "--arch-graphs option is mandatory");

  profile(arch_graphs_stream, options);

  return EXIT_SUCCESS;
}
