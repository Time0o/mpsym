#include <cassert>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

#include "dbg.h"
#include "perm.h"
#include "pr_randomizer.h"
#include "schreier_tree.h"
#include "schreier_sims.h"

namespace cgtl
{

namespace schreier_sims
{

std::vector<std::vector<unsigned>> orbits(std::vector<Perm> const &generators)
{
  unsigned n = generators[0].degree();

  std::vector<std::vector<unsigned>> res;
  std::vector<int> orbit_indices(n + 1u, -1);

  unsigned processed = 0u;

  for (auto i = 1u; i <= n; ++i) {
    int orbit_index1 = orbit_indices[i];
    if (orbit_index1 == -1) {
      orbit_index1 = static_cast<int>(res.size());
      orbit_indices[i] = orbit_index1;

      res.push_back({i});

      if (++processed == n)
        return res;
    }

    for (Perm const &gen : generators) {
      unsigned j = gen[i];

      int orbit_index2 = orbit_indices[j];
      if (orbit_index2 == -1) {
        res[orbit_index1].push_back(j);
        orbit_indices[j] = orbit_index1;

        if (++processed == n)
          return res;
      }
    }
  }

  return res;
}

std::vector<unsigned> orbit(
  unsigned alpha, std::vector<Perm> const &generators, SchreierTree &st)
{
  assert(generators.size() > 0u && "generator set not empty");
  assert(alpha <= generators[0].degree() && "alpha <= N");

  std::vector<unsigned> result {alpha};
  st.create_root(alpha);

  std::vector<unsigned> stack {alpha};
  std::set<unsigned> done {alpha};

  while (!stack.empty()) {
    unsigned beta = stack.back();
    stack.pop_back();

    for (auto i = 0u; i < generators.size(); ++i) {
      Perm const &gen = generators[i];
      unsigned beta_prime = gen[beta];

      if (done.find(beta_prime) == done.end()) {
        result.push_back(beta_prime);
        done.insert(beta_prime);
        stack.push_back(beta_prime);

        st.create_edge(beta_prime, beta, generators[i]);
      }
    }
  }

  return result;
}

std::pair<Perm, unsigned> strip(
  Perm const &perm, std::vector<unsigned> const &base,
  std::vector<SchreierTree> const &sts)
{
  Perm result(perm);

  for (unsigned i = 0u; i < base.size(); ++i) {
    unsigned beta = result[base[i]];
    if (!sts[i].contains(beta))
      return std::make_pair(result, i + 1u);

    result *= ~sts[i].transversal(beta);
  }

  return std::make_pair(result, base.size() + 1u);
}

static void schreier_sims_init(
  std::vector<unsigned> &base, std::vector<Perm> &generators,
  std::vector<SchreierTree> &sts,
  std::vector<std::vector<Perm>> &strong_generators,
  std::vector<std::vector<unsigned>> &fundamental_orbits)
{
  assert(generators.size() > 0u && "generator set not empty");

  unsigned degree = generators[0].degree();

#ifndef NDEBUG
  for (auto const &gen : generators)
    assert(gen.degree() == degree && "all generators have same degree");
#endif

  Dbg(Dbg::DBG) << "=== Input";
  Dbg(Dbg::DBG) << "B = " << base;
  Dbg(Dbg::DBG) << "S = " << generators;

  // add initial base points
  for (auto it  = generators.begin(); it != generators.end(); ) {
    Perm gen = *it;

    if (gen.id()) {
      it = generators.erase(it);

    } else {
      ++it;

      bool stabilizes = true;
      for (unsigned b : base) {
        if (gen[b] != b) {
          stabilizes = false;
          break;
        }
      }

      if (stabilizes) {
#ifndef NDEBUG
        bool extended_base = false;
#endif
        for (unsigned i = 1u; i <= degree; ++i) {
          if (gen[i] != i) {
            base.push_back(i);
#ifndef NDEBUG
            extended_base = true;
#endif
            break;
          }
        }

        assert(extended_base && "no generator fixes all base elements");
      }
    }
  }

  strong_generators.resize(base.size());
  fundamental_orbits.resize(base.size());
  sts.resize(base.size(), SchreierTree(degree));

  // calculate initial strong generator sets
  for (unsigned i = 0u; i < base.size(); ++i) {
    for (Perm const &gen : generators) {
      bool stabilizes = true;
      for (unsigned k = 0u; k < i; ++k) {
        if (gen[base[k]] != base[k]) {
          stabilizes = false;
          break;
        }
      }
      if (stabilizes) {
        strong_generators[i].push_back(gen);
      }
    }

    fundamental_orbits[i] = orbit(base[i], strong_generators[i], sts[i]);
  }

  Dbg(Dbg::DBG) << "=== Initial values";
  Dbg(Dbg::DBG) << "B = " << base;
  for (unsigned i = 0u; i < base.size(); ++i) {
    Dbg(Dbg::DBG) << "S(" << (i + 1u) << ") = " << strong_generators[i];
    Dbg(Dbg::DBG) << "O(" << (i + 1u) << ") = " << fundamental_orbits[i];
  }
}

static void schreier_sims_finish(std::vector<unsigned> const &base,
  std::vector<Perm> &generators,
  std::vector<std::vector<Perm>> const &strong_generators)
{
  std::unordered_set<Perm> unique_generators;

  for (auto const &gens : strong_generators)
    unique_generators.insert(gens.begin(), gens.end());

  generators = std::vector<Perm>(unique_generators.begin(),
                                 unique_generators.end());

  Dbg(Dbg::DBG) << "=== Result";
  Dbg(Dbg::DBG) << "B = " << base;
  Dbg(Dbg::DBG) << "SGS = " << generators;
}

void schreier_sims(std::vector<unsigned> &base,
  std::vector<Perm> &generators, std::vector<SchreierTree> &sts)
{
  Dbg(Dbg::DBG) << "Executing schreier sims algorithm";

  // intialize
  std::vector<std::vector<Perm>> strong_generators;
  std::vector<std::vector<unsigned>> fundamental_orbits;

  schreier_sims_init(
    base, generators, sts, strong_generators, fundamental_orbits);

  unsigned degree = generators[0].degree();

  unsigned i = base.size();
  while (i >= 1u) {
top:
    Dbg(Dbg::TRACE) << "=== Main loop (i = " << i << ")";
    auto const &st = sts[i - 1u];

    for (unsigned beta : fundamental_orbits[i - 1u]) {
      Dbg(Dbg::TRACE) << "== Orbit element " << beta;

      Perm u_beta = st.transversal(beta);

      for (Perm const &x : strong_generators[i - 1u]) {
        unsigned beta_x = x[beta];

        Perm schreier_generator = u_beta * x * ~st.transversal(beta_x);
        Dbg(Dbg::TRACE) << "Schreier Generator: " << schreier_generator;

        if (schreier_generator.id())
          continue;

        bool update_strong_generators = false;
        std::pair<Perm, unsigned> strip_result =
          strip(schreier_generator, base, sts);

        Perm strip_perm = std::get<0>(strip_result);
        unsigned strip_level = std::get<1>(strip_result);
        Dbg(Dbg::TRACE) << "Strips to: " << strip_perm << ", " << strip_level;

        if (strip_level <= base.size()) {
           update_strong_generators = true;

        } else if (strip_perm != Perm(degree)) {
          update_strong_generators = true;

          unsigned next_basepoint = 1u;
          while (strip_perm[next_basepoint] == next_basepoint)
            ++next_basepoint;

          base.push_back(next_basepoint);

          Dbg(Dbg::TRACE) << "Adjoined new basepoint:";
          Dbg(Dbg::TRACE) << "B = " << base;
        }

        if (update_strong_generators) {
          Dbg(Dbg::TRACE) << "Updating strong generators:";

          for (unsigned j = i; j < strip_level; ++ j) {
            if (strong_generators.size() <= j) {
              strong_generators.push_back({});
              fundamental_orbits.push_back({});
              sts.push_back(SchreierTree(degree));
            }

            strong_generators[j].push_back(strip_perm);

            fundamental_orbits[j] =
              orbit(base[j], strong_generators[j], sts[j]);

            Dbg(Dbg::TRACE) << "S(" << (j + 1u) << ") = " << strong_generators[j];
            Dbg(Dbg::TRACE) << "O(" << (j + 1u) << ") = " << fundamental_orbits[j];
          }

          i = strip_level;
          goto top;
        }
      }
    }
    i -= 1u;
  }

  schreier_sims_finish(base, generators, strong_generators);
}

void schreier_sims_random(
  std::vector<unsigned> &base, std::vector<Perm> &generators,
  std::vector<SchreierTree> &sts, unsigned w)
{
  Dbg(Dbg::DBG) << "Executing (random) schreier sims algorithm";

  // intialize
  std::vector<std::vector<Perm>> strong_generators;
  std::vector<std::vector<unsigned>> fundamental_orbits;

  schreier_sims_init(
    base, generators, sts, strong_generators, fundamental_orbits);

  unsigned degree = generators[0].degree();

  // random group element generator
  PrRandomizer pr(generators);

  unsigned c = 0u;
  while (c < w) {
    Perm rand_perm = pr.next();
    Dbg(Dbg::TRACE) << "Random group element: " << rand_perm;

    bool update_strong_generators = false;
    std::pair<Perm, unsigned> strip_result = strip(rand_perm, base, sts);

    Perm strip_perm = std::get<0>(strip_result);
    unsigned strip_level = std::get<1>(strip_result);
    Dbg(Dbg::TRACE) << "Strips to: " << strip_perm << ", " << strip_level;

    if (strip_level <= base.size()) {
      update_strong_generators = true;

    } else if (!strip_perm.id()) {
      update_strong_generators = true;

      for (unsigned i = 1u; i <= degree; ++i) {
        if (strip_perm[i] != i) {
          base.push_back(i);
          strong_generators.push_back({});

          Dbg(Dbg::TRACE) << "Adjoined new basepoint:";
          Dbg(Dbg::TRACE) << "B = " << base;

          break;
        }
      }
    }

    if (update_strong_generators) {
      Dbg(Dbg::TRACE) << "Updating strong generators:";

      strong_generators.resize(strip_level);
      fundamental_orbits.resize(strip_level);
      sts.resize(strip_level, SchreierTree(degree));

      for (unsigned i = 1u; i < strip_level; ++i) {
        strong_generators[i].push_back(strip_perm);
        fundamental_orbits[i] = orbit(base[i], strong_generators[i], sts[i]);

        Dbg(Dbg::TRACE) << "S(" << (i + 1u) << ") = " << strong_generators[i];
        Dbg(Dbg::TRACE) << "O(" << (i + 1u) << ") = " << fundamental_orbits[i];
      }

      c = 0u;

    } else { ++c; }
  }

  schreier_sims_finish(base, generators, strong_generators);
}

} // namespace schreier_sims

} // namespace cgtl
