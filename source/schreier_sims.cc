#include <cassert>
#include <memory>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

#include "dbg.h"
#include "perm.h"
#include "pr_randomizer.h"
#include "schreier_sims.h"
#include "schreier_structure.h"

/**
 * @file schreier_sims.cc
 * @brief Implements data structures and functions defined in schreier_sims.h.
 *
 * @author Timo Nicolai
 */

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
  unsigned alpha, std::vector<Perm> const &generators,
  std::shared_ptr<SchreierStructure> st)
{
  assert(generators.size() > 0u && "generator set not empty");
  assert(alpha <= generators[0].degree() && "alpha <= N");

  std::vector<unsigned> result {alpha};

  if (st) {
    st->create_root(alpha);
    st->create_labels(generators);
  }

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

        if (st)
          st->create_edge(beta_prime, beta, i);
      }
    }
  }

  return result;
}

std::pair<Perm, unsigned> strip(
  Perm const &perm, std::vector<unsigned> const &base,
  std::vector<std::shared_ptr<SchreierStructure>> const &sts)
{
  Perm result(perm);

  for (unsigned i = 0u; i < base.size(); ++i) {
    unsigned beta = result[base[i]];
    if (!sts[i]->contains(beta))
      return std::make_pair(result, i + 1u);

    result *= ~sts[i]->transversal(beta);
  }

  return std::make_pair(result, base.size() + 1u);
}

void schreier_sims_finish(
  std::vector<unsigned> const &base, std::vector<Perm> &generators,
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

} // namespace schreier_sims

} // namespace cgtl
