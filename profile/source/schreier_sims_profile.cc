#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#include <getopt.h>
#include <libgen.h>

#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"
#include "permlib.h"
#include "timer.h"
#include "util.h"

#include "profile_parse.h"
#include "profile_args.h"
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
    "-i|--implementation  {gap|mpsym|permlib}",
    "[-s|--schreier-sims] {deterministic|random|random-no-guarantee}",
    "[-t|--transversals]  {explicit|schreier-trees|shallow-schreier-trees}",
    "[--check-altsym]",
    "[--no-reduce-gens]",
    "[-g|--groups GROUPS]",
    "[-a|--arch-graph ARCH_GRAPH]",
    "[-r|--num-runs]",
    "[--num-discarded-runs NUM_DISCARDED_RUNS]",
    "[--summarize-runs]",
    "[-v|--verbose]",
    "[--show-gap-errors]"
  };

  s << "usage: " << progname << '\n';
  for (char const *opt : opts)
    s << "  " << opt << '\n';
}

struct ProfileOptions
{
  VariantOption implementation{"gap", "mpsym", "permlib"};
  VariantOption schreier_sims{"deterministic", "random", "random-no-guarantee"};
  VariantOption transversals{"explicit", "schreier-trees", "shallow-schreier-trees"};
  bool check_altsym = false;
  bool no_reduce_gens = false;
  bool groups_input = false;
  bool arch_graph_input = false;
  unsigned num_runs = 1u;
  unsigned num_discarded_runs = 0u;
  bool summarize_runs = false;
  bool verbose = false;
  bool show_gap_errors = false;
};

cgtl::BSGS::Options bsgs_options_mpsym(ProfileOptions const &options)
{
  using cgtl::BSGS;

  BSGS::Options bsgs_options;

  if (options.schreier_sims.is("deterministic")) {
    bsgs_options.construction = BSGS::Construction::SCHREIER_SIMS;
  } else if (options.schreier_sims.is("random")) {
    bsgs_options.construction = BSGS::Construction::SCHREIER_SIMS_RANDOM;
  } else if (options.schreier_sims.is("random-no-guarantee")) {
    bsgs_options.construction = BSGS::Construction::SCHREIER_SIMS_RANDOM;
    bsgs_options.schreier_sims_random_guarantee = false;
  } else {
    throw std::logic_error("unreachable");
  }

  if (options.transversals.is("explicit"))
    bsgs_options.transversals = BSGS::Transversals::EXPLICIT;
  else if (options.transversals.is("schreier-trees"))
    bsgs_options.transversals = BSGS::Transversals::SCHREIER_TREES;
  else if (options.transversals.is("shallow-schreier-trees"))
    bsgs_options.transversals = BSGS::Transversals::SHALLOW_SCHREIER_TREES;
  else
    throw std::logic_error("unreachable");

  if (options.check_altsym)
    bsgs_options.check_altsym = true;

  if (options.no_reduce_gens)
    bsgs_options.reduce_gens = false;

  return bsgs_options;
}

std::string make_perm_group_gap(gap::PermSet const &generators,
                                ProfileOptions const &options)
{
  (void)options;

  return "StabChain(Group(" + generators.permutations + "));\n";
}

void make_perm_group_mpsym(cgtl::PermSet const &generators,
                           ProfileOptions const &options)
{
  using cgtl::BSGS;
  using cgtl::PermGroup;
  using cgtl::PermSet;

  auto bsgs_options(bsgs_options_mpsym(options));

  PermGroup g(BSGS(generators.degree(), generators, &bsgs_options));
}

template <typename T>
struct TypeTag { using type = T; };

std::variant<
  TypeTag<permlib::ExplicitTransversal<permlib::Permutation>>,
  TypeTag<permlib::SchreierTreeTransversal<permlib::Permutation>>,
  TypeTag<permlib::ShallowSchreierTreeTransversal<permlib::Permutation>>
>
permlib_transversal_type(VariantOption const &transversals)
{
  using namespace permlib;

  if (transversals.is("explicit"))
    return TypeTag<ExplicitTransversal<Permutation>>{};
  else if (transversals.is("schreier-trees"))
    return TypeTag<SchreierTreeTransversal<Permutation>>{};
  else if (transversals.is("shallow-schreier-trees"))
    return TypeTag<ShallowSchreierTreeTransversal<Permutation>>{};
  else
    throw std::logic_error("unreachable");
}

void make_perm_group_permlib(permlib::PermSet const &generators,
                             ProfileOptions const &options)
{
  using namespace permlib;

  std::visit([&](auto transv_type_) {
    using transv_type = typename decltype(transv_type_)::type;

    if (options.schreier_sims.is("deterministic")) {
      SchreierSimsConstruction<Permutation, transv_type>
      construction(generators.degree);

      construction.construct(generators.permutations.begin(),
                             generators.permutations.end());

    } else if (options.schreier_sims.is("random")) {
      BSGS<Permutation, transv_type> bsgs(generators.degree);

      std::unique_ptr<BSGSRandomGenerator<Permutation, transv_type>>
      random_generator(new BSGSRandomGenerator<Permutation, transv_type>(bsgs));

      RandomSchreierSimsConstruction<Permutation, transv_type>
      construction(generators.degree, random_generator.get());

      bool guaranteed = true;
      construction.construct(generators.permutations.begin(),
                             generators.permutations.end(),
                             guaranteed);
    }
  }, permlib_transversal_type(options.transversals));
}

std::vector<double> run_groups(unsigned degree,
                               std::string const &generators,
                               ProfileOptions const &options)
{
  std::vector<double> ts;

  if (options.implementation.is("gap")) {
    auto generators_gap(parse_generators_gap(degree, generators));

    run_gap({},
            make_perm_group_gap(generators_gap, options),
            options.num_discarded_runs,
            options.num_runs,
            options.verbose,
            !options.show_gap_errors,
            &ts);

  } else if (options.implementation.is("mpsym")) {
    auto generators_mpsym(parse_generators_mpsym(degree, generators));

    run_cpp([&]{ make_perm_group_mpsym(generators_mpsym, options); },
            options.num_discarded_runs,
            options.num_runs,
            &ts);

  } else if (options.implementation.is("permlib")) {
    auto generators_permlib(parse_generators_permlib(degree, generators));

    run_cpp([&]{ make_perm_group_permlib(generators_permlib, options); },
            options.num_discarded_runs,
            options.num_runs,
            &ts);
  }

  return ts;
}

std::vector<double> run_arch_graph(std::string const &arch_graph,
                                   ProfileOptions const &options)
{
  using boost::multiprecision::cpp_int;

  using cgtl::ArchGraphSystem;
  using cgtl::PermSet;

  std::vector<double> ts;

  auto ag(ArchGraphSystem::from_lua(arch_graph));

  cpp_int num_automorphisms;
  PermSet automorphisms_generators;

  if (options.implementation.is("gap")) {
    auto gap_script("G:=" + ag->to_gap() + ";\n"
                    "StabChain(G);\n"
                    "Print(Size(G), \";\\n\");\n"
                    "Print(GeneratorsOfGroup(G), \";\\n\");\n");

    auto gap_output(run_gap({"grape"},
                            gap_script,
                            options.num_discarded_runs,
                            options.num_runs,
                            true,
                            !options.show_gap_errors,
                            &ts));

    num_automorphisms = cpp_int(gap_output[0]);
    automorphisms_generators = parse_generators_mpsym(0, gap_output[1]);

  } else if (options.implementation.is("mpsym")) {
    auto bsgs_options(bsgs_options_mpsym(options));

    run_cpp([&]{ ag->invalidate_automorphisms();
                 ag->automorphisms(&bsgs_options); },
            options.num_discarded_runs,
            options.num_runs,
            &ts);

    num_automorphisms = ag->automorphisms().order();
    automorphisms_generators = ag->automorphisms_generators();

  } else if (options.implementation.is("permlib")) {
    throw std::logic_error("graph automorphisms not supported by permlib");
  }

  if (options.verbose) {
    info("Automorphism group has order:");
    info(num_automorphisms);

    info("Automorphism generators are:");
    info(automorphisms_generators);
  }

  return ts;
}

template<typename FUNC>
void run(FUNC &&f,
         ProfileOptions const &options)
{
  auto ts(f(options));

  if (options.summarize_runs) {
    double t_mean, t_stddev;
    util::mean_stddev(ts, &t_mean, &t_stddev);

    result("Mean:", t_mean, "s");
    result("Stddev:", t_stddev, "s");

  } else {
    result("Runtimes:");
    for (double t : ts)
      result(t, "s");
  }

  if (options.verbose && options.implementation.is("mpsym")) {
    debug("Timer dumps:");
    debug_timer_dump("strip");
    debug_timer_dump("extend base");
    debug_timer_dump("update strong gens");
  }
}

void do_profile(Stream &instream,
                ProfileOptions const &options)
{
  using namespace std::placeholders;

  if (options.verbose) {
    debug("Implementation:", options.implementation.get());

    if (!options.implementation.is("gap")) {
      debug("Schreier-sims variant:", options.schreier_sims.get());
      debug("Transversals:", options.transversals.get());
    }
  }

  if (options.groups_input) {
    foreach_line(instream.stream, [&](std::string const &line, unsigned lineno){
      auto group(parse_group(line));

      if (options.verbose) {
        info("Constructing group", lineno);
        info("=> degree", group.degree);
        info("=> orders", group.order);
        info("=> generators", group.generators);
      } else {
        info("Constructing group", lineno);
      }

      run(std::bind(run_groups, group.degree, group.generators, _1), options);
    });

  } else if (options.arch_graph_input) {
    if (options.verbose)
      info("Constructing automorphism group");

    run(std::bind(run_arch_graph, read_file(instream.stream), _1), options);
  }
}

} // anonymous namespace

int main(int argc, char **argv)
{
  progname = basename(argv[0]);

  struct option long_options[] = {
    {"help",                no_argument,       0,       'h'},
    {"implementation",      required_argument, 0,       'i'},
    {"schreier-sims",       required_argument, 0,       's'},
    {"transversals",        required_argument, 0,       't'},
    {"check-altsym",        no_argument,       0,        1 },
    {"no-reduce-gens",      no_argument,       0,        2 },
    {"groups",              no_argument,       0,       'g'},
    {"arch-graph",          no_argument,       0,       'a'},
    {"num-runs",            required_argument, 0,       'r'},
    {"num-discarded-runs",  required_argument, 0,        3 },
    {"summarize-runs",      no_argument, 0,              4 },
    {"verbose",             no_argument,       0,       'v'},
    {"show-gap-errors",     no_argument,       0,        5 },
    {nullptr,               0,                 nullptr,  0 }
  };

  ProfileOptions options;

  Stream instream;

  for (;;) {
    int c = getopt_long(argc, argv, "hi:s:t:g:a:r:v", long_options, nullptr);
    if (c == -1)
      break;

    try {
      switch(c) {
      case 'h':
        usage(std::cout);
        return EXIT_SUCCESS;
      case 'i':
        options.implementation.set(optarg);
        break;
      case 's':
        options.schreier_sims.set(optarg);
        break;
      case 't':
        options.transversals.set(optarg);
        break;
      case 1:
        options.check_altsym = true;
        break;
      case 2:
        options.no_reduce_gens = true;
        break;
      case 'g':
        OPEN_STREAM(instream, optarg);
        options.groups_input = true;
        break;
      case 'a':
        OPEN_STREAM(instream, optarg);
        options.arch_graph_input = true;
        break;
      case 'r':
        options.num_runs = stox<unsigned>(optarg);
        break;
      case 3:
        options.num_discarded_runs = stox<unsigned>(optarg);
        break;
      case 4:
        options.summarize_runs = true;
        break;
      case 'v':
        options.verbose = true;
        TIMER_ENABLE();
        break;
      case 5:
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

  CHECK_OPTION(options.implementation.is_set(),
               "--implementation option is mandatory");

  CHECK_OPTION((options.implementation.is("gap") || options.schreier_sims.is_set()),
               "--schreier-sims option is mandatory when not using gap");

  CHECK_OPTION((options.implementation.is("gap") || options.transversals.is_set()),
               "--transversal-storage option is mandatory when not using gap");

  CHECK_OPTION(options.groups_input != options.arch_graph_input,
               "EITHER --arch-graph OR --groups must be given");

  try {
    do_profile(instream, options);
  } catch (std::exception const &e) {
    error("profiling failed:", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
