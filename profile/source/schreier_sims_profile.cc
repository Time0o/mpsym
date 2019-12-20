#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

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

namespace
{

std::string progname;

void usage(std::ostream &s)
{
  char const *opts[] = {
    "[-h|--help]",
    "-i|--implementation        {gap|mpsym|permlib}",
    "[-s|--schreier-sims]       {deterministic|random}",
    "[-t|--transversal-storage] {explicit|schreier-trees|shallow-schreier-trees}",
    "[-c|--num-cycles]",
    "[-r|--num-runs]",
    "[--realtime-clock]",
    "[-v|--verbose]",
    "GROUPS"
  };

  s << "usage: " << progname << '\n';
  for (char const *opt : opts)
    s << "  " << opt << '\n';
}

struct ProfileOptions
{
  VariantOption library{"gap", "mpsym", "permlib"};
  VariantOption schreier_sims{"deterministic", "random"};
  VariantOption transversals{"explicit", "schreier-trees", "shallow-schreier-trees"};
  unsigned num_cycles = 1u;
  unsigned num_runs = 1u;
  bool verbose = false;
};

std::string make_perm_group_gap(gap::PermSet const &generators,
                                ProfileOptions const &options)
{
  std::stringstream ss;
  ss << "for i in [1.." << options.num_cycles << "] do\n";
  ss << "  StabChain(Group(" + generators.permutations + "));\n";
  ss << "od;\n";

  return ss.str();
}

void make_perm_group_mpsym(cgtl::PermSet const &generators,
                           ProfileOptions const &options)
{
  using cgtl::BSGS;
  using cgtl::PermGroup;
  using cgtl::PermSet;

  BSGS::Construction constr;
  if (options.schreier_sims.is("deterministic"))
    constr = BSGS::Construction::SCHREIER_SIMS;
  else if (options.schreier_sims.is("random"))
    constr = BSGS::Construction::SCHREIER_SIMS_RANDOM;
  else
    throw std::logic_error("unreachable");

  BSGS::Transversals transv;
  if (options.transversals.is("explicit"))
    transv = BSGS::Transversals::EXPLICIT;
  else if (options.transversals.is("schreier-trees"))
    transv = BSGS::Transversals::SCHREIER_TREES;
  else if (options.transversals.is("shallow-schreier-trees"))
    transv = BSGS::Transversals::SHALLOW_SCHREIER_TREES;
  else
    throw std::logic_error("unreachable");

  for (unsigned i = 0u; i < options.num_cycles; ++i)
    PermGroup g(generators.degree(), generators, constr, transv);
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

      for (unsigned i = 0; i < options.num_cycles; ++i)
        construction.construct(generators.permutations.begin(),
                               generators.permutations.end());

    } else if (options.schreier_sims.is("random")) {
      BSGS<Permutation, transv_type> bsgs(generators.degree);

      std::unique_ptr<BSGSRandomGenerator<Permutation, transv_type>>
      random_generator(new BSGSRandomGenerator<Permutation, transv_type>(bsgs));

      RandomSchreierSimsConstruction<Permutation, transv_type>
      construction(generators.degree, random_generator.get());

      bool guaranteed = true;
      for (unsigned i = 0; i < options.num_cycles; ++i)
        construction.construct(generators.permutations.begin(),
                               generators.permutations.end(),
                               guaranteed);
    }
  }, permlib_transversal_type(options.transversals));
}

std::vector<double> run(std::string const &generators,
                        ProfileOptions const &options)
{
  std::vector<double> ts;
  for (unsigned r = 0; r < options.num_runs; ++r) {
    if (options.verbose)
      debug_progress("Executing run", r + 1, "/", options.num_runs);

    double t;
    if (options.library.is("gap")) {
      run_gap(make_perm_group_gap(parse_generators_gap(generators), options),
              options.verbose,
              &t);
    } else {
      if (options.library.is("mpsym")) {
        run_cpp([&]{
          make_perm_group_mpsym(parse_generators_mpsym(generators), options);
        }, &t);
      } else if (options.library.is("permlib")) {
        run_cpp([&]{
          make_perm_group_permlib(parse_generators_permlib(generators), options);
        }, &t);
      }
    }

    ts.push_back(t);
  }

  if (options.verbose) {
    debug_progress_done();

    if (options.library.is("mpsym")) {
      debug("Timer dumps:");
      debug_timer_dump("strip");
      debug_timer_dump("extend base");
      debug_timer_dump("update strong gens");
    }
  }

  return ts;
}

void profile(std::ifstream &groups_stream,
             ProfileOptions const &options)
{
  if (options.verbose) {
    debug("Implementation:", options.library.get());
    debug("Schreier-sims variant:", options.schreier_sims.get());
    debug("Transversals:", options.transversals.get());

    if (options.num_cycles > 1)
      debug("Constructions per run:", options.num_cycles);
  }

  foreach_line(groups_stream, [&](std::string const &line, unsigned lineno){
    auto group(parse_group(line));

    if (options.verbose) {
      info("Constructing group", lineno,
           "with degree", group.degree,
           "and generators", group.generators);
    } else {
      info("Constructing group", lineno,
           "with degree", group.degree);
    }

    auto ts = run(group.generators, options);

    double t_mean, t_stddev;
    util::mean_stddev(ts, &t_mean, &t_stddev);

    result("Mean:", t_mean, "s");
    result("Stddev:", t_stddev, "s");
  });
}

} // namespace

int main(int argc, char **argv)
{
  progname = basename(argv[0]);

  struct option long_options[] = {
    {"help",                no_argument,       0,       'h'},
    {"implementation",      required_argument, 0,       'i'},
    {"schreier-sims",       required_argument, 0,       's'},
    {"transversal-storage", required_argument, 0,       't'},
    {"num-cyles",           required_argument, 0,       'c'},
    {"num-runs",            required_argument, 0,       'r'},
    {"realtime-clock",      no_argument,       0,        1 },
    {"verbose",             no_argument,       0,       'v'},
    {nullptr,               0,                 nullptr,  0 }
  };

  ProfileOptions options;

  std::ifstream groups_stream;

  for (;;) {
    int c = getopt_long(argc, argv, "hi:s:t:r:c:v", long_options, nullptr);
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
      case 's':
        options.schreier_sims.set(optarg);
        break;
      case 't':
        options.transversals.set(optarg);
        break;
      case 'c':
        options.num_cycles = stox<unsigned>(optarg);
        break;
      case 'r':
        options.num_runs = stox<unsigned>(optarg);
        break;
      case 1:
        timer_realtime_enable();
        break;
      case 'v':
        options.verbose = true;
        TIMER_ENABLE();
        break;
      default:
        return EXIT_FAILURE;
      }
    } catch (std::invalid_argument const &) {
      error("invalid option argument");
      return EXIT_FAILURE;
    }
  }

  CHECK_OPTION(options.library.is_set(),
               "--implementation option is mandatory");

  CHECK_OPTION((options.library.is("gap") || options.schreier_sims.is_set()),
               "--schreier-sims option is mandatory when not using gap");

  CHECK_OPTION((options.library.is("gap") || options.transversals.is_set()),
               "--transversal-storage option is mandatory when not using gap");

  CHECK_ARGUMENT("GROUPS");

  OPEN_STREAM(groups_stream, argv[optind]);

  try {
    profile(groups_stream, options);
  } catch (std::exception const &e) {
    error("profiling failed:", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
