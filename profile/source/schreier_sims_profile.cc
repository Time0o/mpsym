#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
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
    "-i|--implementation        {gap|mpsym|permlib}",
    "[-s|--schreier-sims]       {deterministic|random}",
    "[-t|--transversal-storage] {explicit|schreier-trees|shallow-schreier-trees}",
    "[-c|--num-cycles]",
    "[-r|--num-runs]",
    "[-v|--verbose]",
    "GROUPS"
  };

  s << "usage: " << progname << '\n';
  for (char const *opt : opts)
    s << "  " << opt << '\n';
}

void make_perm_group_mpsym(VariantOption const &schreier_sims_impl,
                           VariantOption const &transversals_impl,
                           unsigned degree,
                           cgtl::PermSet const &gens,
                           unsigned num_cycles)
{
  using cgtl::BSGS;
  using cgtl::PermGroup;
  using cgtl::PermSet;

  BSGS::Construction constr;
  if (schreier_sims_impl.is("deterministic"))
    constr = BSGS::Construction::SCHREIER_SIMS;
  else if (schreier_sims_impl.is("random"))
    constr = BSGS::Construction::SCHREIER_SIMS_RANDOM;
  else
    throw std::logic_error("unreachable");

  BSGS::Transversals transv;
  if (transversals_impl.is("explicit"))
    transv = BSGS::Transversals::EXPLICIT;
  else if (transversals_impl.is("schreier-trees"))
    transv = BSGS::Transversals::SCHREIER_TREES;
  else if (transversals_impl.is("shallow-schreier-trees"))
    transv = BSGS::Transversals::SHALLOW_SCHREIER_TREES;
  else
    throw std::logic_error("unreachable");

  for (unsigned i = 0u; i < num_cycles; ++i)
    PermGroup g(degree, gens, constr, transv);
}

template <typename T>
struct TypeTag { using type = T; };

std::variant<
  TypeTag<permlib::ExplicitTransversal<permlib::Permutation>>,
  TypeTag<permlib::SchreierTreeTransversal<permlib::Permutation>>,
  TypeTag<permlib::ShallowSchreierTreeTransversal<permlib::Permutation>>
>
permlib_transversal_type(VariantOption const &transversals_impl)
{
  using namespace permlib;

  if (transversals_impl.is("explicit"))
    return TypeTag<ExplicitTransversal<Permutation>>{};
  else if (transversals_impl.is("schreier-trees"))
    return TypeTag<SchreierTreeTransversal<Permutation>>{};
  else if (transversals_impl.is("shallow-schreier-trees"))
    return TypeTag<ShallowSchreierTreeTransversal<Permutation>>{};
  else
    throw std::logic_error("unreachable");
}

void make_perm_group_permlib(VariantOption const &schreier_sims_impl,
                             VariantOption const &transversals_impl,
                             unsigned degree,
                             std::vector<permlib::Permutation::ptr> const &gens,
                             unsigned num_cycles)
{
  using namespace permlib;

  std::visit([&](auto transv_type_) {
    using transv_type = typename decltype(transv_type_)::type;

    if (schreier_sims_impl.is("deterministic")) {
      SchreierSimsConstruction<Permutation, transv_type>
      construction(degree);

      for (unsigned i = 0; i < num_cycles; ++i)
        construction.construct(gens.begin(), gens.end());

    } else if (schreier_sims_impl.is("random")) {
      BSGS<Permutation, transv_type> bsgs(degree);

      std::unique_ptr<BSGSRandomGenerator<Permutation, transv_type>>
      random_generator(new BSGSRandomGenerator<Permutation, transv_type>(bsgs));

      RandomSchreierSimsConstruction<Permutation, transv_type>
      construction(degree, random_generator.get());

      bool guaranteed = true;
      for (unsigned i = 0; i < num_cycles; ++i)
        construction.construct(gens.begin(), gens.end(), guaranteed);
    }
  }, permlib_transversal_type(transversals_impl));
}

std::vector<double> run(VariantOption const &library_impl,
                        VariantOption const &schreier_sims_impl,
                        VariantOption const &transversals_impl,
                        unsigned degree,
                        std::string const &generators_str,
                        unsigned num_cycles,
                        unsigned num_runs,
                        bool verbose)
{
  std::vector<double> ts;
  for (unsigned r = 0; r < num_runs; ++r) {
    if (verbose)
      progress("Executing run", r + 1, "/", num_runs);

    double t;
    if (library_impl.is("gap")) {
      auto generators(parse_generators_gap(generators_str));

      std::stringstream ss;
      ss << "for i in [1.." << num_cycles << "] do\n";
      ss << "  StabChain(Group(" + generators + "));\n";
      ss << "od;\n";

      run_gap(ss.str(), &t);

    } else {
      if (library_impl.is("mpsym")) {
        auto generators(parse_generators_mpsym(generators_str));

        run_cpp([&]{
          make_perm_group_mpsym(schreier_sims_impl,
                                transversals_impl,
                                degree,
                                generators,
                                num_cycles);
          }, &t);

      } else if (library_impl.is("permlib")) {
        auto generators(parse_generators_permlib(generators_str));

        run_cpp([&]{
          make_perm_group_permlib(schreier_sims_impl,
                                  transversals_impl,
                                  degree,
                                  generators,
                                  num_cycles);
          }, &t);
      }
    }

    ts.push_back(t);
  }

  if (verbose)
    progress_done();

  return ts;
}

void profile(VariantOption const &library_impl,
             VariantOption const &schreier_sims_impl,
             VariantOption const &transversals_impl,
             std::ifstream &groups_stream,
             unsigned num_cycles,
             unsigned num_runs,
             bool verbose)
{
  if (verbose) {
    Timer::enabled = true;

    info("Implementation:", library_impl.get());
    info("Schreier-sims variant:", schreier_sims_impl.get());
    info("Transversals:", transversals_impl.get());

    if (num_cycles > 1)
      info("Constructions per run:", num_cycles);
  }

  std::string line;
  int lineno = 1;
  while (std::getline(groups_stream, line)) {
    unsigned degree;
    unsigned order;
    std::string generators_str;

    std::tie(degree, order, generators_str) = parse_group(line);

    if (verbose) {
      info("Constructing group", lineno,
           "with degree", degree,
           "and order", order);
    }

    auto ts = run(library_impl,
                  schreier_sims_impl,
                  transversals_impl,
                  degree,
                  generators_str,
                  num_cycles,
                  num_runs,
                  verbose);

    double t_mean, t_stddev;
    util::mean_stddev(ts, &t_mean, &t_stddev);

    result("Mean:", t_mean, "s");
    result("Stddev:", t_stddev, "s");

    ++lineno;
  }

  if (groups_stream.bad())
    throw std::runtime_error("failed to read groups file");
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

  VariantOption library_impl({"gap", "mpsym", "permlib"});

  VariantOption schreier_sims_impl({"deterministic", "random"});

  VariantOption transversals_impl(
    {"explicit", "schreier-trees", "shallow-schreier-trees"});

  unsigned num_cycles = 1;
  unsigned num_runs = 1;
  bool verbose = false;

  // TODO: detect superfluous non-argument options

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
        library_impl.set(optarg);
        break;
      case 's':
        schreier_sims_impl.set(optarg);
        break;
      case 't':
        transversals_impl.set(optarg);
        break;
      case 'c':
        num_cycles = stox<unsigned>(optarg);
        break;
      case 'r':
        num_runs = stox<unsigned>(optarg);
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

  CHECK_OPTION((library_impl.is("gap") || schreier_sims_impl.is_set()),
               "--schreier-sims option is mandatory when not using gap");

  CHECK_OPTION((library_impl.is("gap") || transversals_impl.is_set()),
               "--transversal-storage option is mandatory when not using gap");

  CHECK_FILE_ARGUMENT("GROUPS");

  std::ifstream groups_stream(argv[optind]);

  try {
    profile(library_impl,
            schreier_sims_impl,
            transversals_impl,
            groups_stream,
            num_cycles,
            num_runs,
            verbose);
  } catch (std::exception const &e) {
    error("profiling failed", e.what());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
