#include <cassert>
#include <set>
#include <utility>
#include <vector>

#include "dbg.h"
#include "perm.h"
#include "schreier_tree.h"
#include "schreier_sims.h"

namespace cgtl
{

std::vector<unsigned> SchreierSims::orbit(unsigned alpha,
  std::vector<Perm> const &generators, SchreierTree &st)
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

std::pair<Perm, unsigned> SchreierSims::strip(Perm const &perm,
  std::vector<unsigned> const &base, std::vector<SchreierTree> const &sts)
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

void SchreierSims::schreier_sims(std::vector<unsigned> &base,
  std::vector<Perm> &generators, std::vector<SchreierTree> &sts)
{
  assert(generators.size() > 0u && "generator set not empty");

  Dbg(Dbg::DBG) << "Executing schreier sims algorithm";
  Dbg(Dbg::DBG) << "=== Input";
  Dbg(Dbg::DBG) << "B = " << base;
  Dbg(Dbg::DBG) << "S = " << generators;

  unsigned degree = generators[0].degree();

  // find initial base point if given base is empty
  if (base.size() == 0u) {
    bool init = false;
    for (unsigned beta = 1u; beta <= degree; ++beta) {
      for (auto const &gen : generators) {
        if (gen[beta] != beta) {
          base.push_back(beta);

          init = true;
          break;
        }
      }
      if (init)
        break;
    }
  }

  std::vector<std::vector<Perm>> strong_generators(base.size());
  std::vector<std::vector<unsigned>> fundamental_orbits(base.size());
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

    fundamental_orbits[i] = orbit(
      base[i], strong_generators[i], sts[i]);
  }

  Dbg(Dbg::DBG) << "=== Initial values";
  Dbg(Dbg::DBG) << "B = " << base;
  for (unsigned i = 0u; i < base.size(); ++i) {
    Dbg(Dbg::DBG) << "S(" << (i + 1u) << ") = " << strong_generators[i];
    Dbg(Dbg::DBG) << "O(" << (i + 1u) << ") = " << fundamental_orbits[i];
  }

  unsigned i = base.size();
  while (i >= 1u) {
top:
    Dbg(Dbg::DBG) << "=== Main loop (i = " << i << ")";
    auto const &st = sts[i - 1u];

    for (unsigned beta : fundamental_orbits[i - 1u]) {
      Dbg(Dbg::DBG) << "== Orbit element " << beta;

      Perm u_beta = st.transversal(beta);

      for (Perm const &x : strong_generators[i - 1u]) {
        unsigned beta_x = x[beta];

        Perm schreier_generator = u_beta * x * ~st.transversal(beta_x);
        Dbg(Dbg::DBG) << "Schreier Generator: " << schreier_generator;

        if (schreier_generator.id())
          continue;

        bool update_strong_generators = false;
        std::pair<Perm, unsigned> strip_result =
          strip(schreier_generator, base, sts);

        Perm strip_perm = std::get<0>(strip_result);
        unsigned strip_level = std::get<1>(strip_result);
        Dbg(Dbg::DBG) << "Strips to: " << strip_perm << ", " << strip_level;

        if (strip_level <= base.size()) {
           update_strong_generators = true;

        } else if (strip_perm != Perm(degree)) {
          update_strong_generators = true;

          unsigned next_basepoint = 1u;
          while (strip_perm[next_basepoint] == next_basepoint)
            ++next_basepoint;

          base.push_back(next_basepoint);

          Dbg(Dbg::DBG) << "Adjoined new basepoint:";
          Dbg(Dbg::DBG) << "B = " << base;
        }

        if (update_strong_generators) {
          Dbg(Dbg::DBG) << "Updating strong generators:";

          for (unsigned j = i; j < strip_level; ++ j) {
            if (strong_generators.size() <= j) {
              strong_generators.push_back({});
              fundamental_orbits.push_back({});
              sts.push_back(SchreierTree(degree));
            }

            strong_generators[j].push_back(strip_perm);

            fundamental_orbits[j] =
              orbit(base[j], strong_generators[j], sts[j]);

            Dbg(Dbg::DBG) << "S(" << (j + 1u) << ") = " << strong_generators[j];
            Dbg(Dbg::DBG) << "O(" << (j + 1u) << ") = " << fundamental_orbits[j];
          }

          i = strip_level;
          goto top;
        }
      }
    }
    i -= 1u;
  }

  generators.clear();
  for (auto const &gens : strong_generators)
    generators.insert(generators.end(), gens.begin(), gens.end());
}

} // namespace cgtl
