#include <algorithm>
#include <cassert>
#include <memory>
#include <ostream>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "perm.h"
#include "schreier_structure.h"

namespace cgtl
{

BSGS::BSGS(unsigned degree,
           PermSet const &generators,
           Construction construction,
           Transversals transversals)
: _degree(degree),
  _transversals(transversals)
{
  generators.assert_degree(degree);

  if (generators.empty() || (generators.size() == 1u && generators[0].id()))
    return;

  switch (construction) {
    case CONSTRUCTION_SCHREIER_SIMS:
      schreier_sims(generators);
      break;
    case CONSTRUCTION_SCHREIER_SIMS_RANDOM:
      schreier_sims_random(generators);
      break;
    case CONSTRUCTION_SOLVE:
      solve(generators);
      break;
    case CONSTRUCTION_AUTO:
      schreier_sims(generators); // TODO
      break;
  }
}

std::vector<unsigned> BSGS::orbit(unsigned i) const
{
  return schreier_structures[i]->nodes();
}

Perm BSGS::transversal(unsigned i, unsigned o) const
{
  return schreier_structures[i]->transversal(o);
}

PermSet BSGS::transversals(unsigned i) const
{
  PermSet transversals;
  for (unsigned o : orbit(i))
    transversals.insert(schreier_structures[i]->transversal(o));

  return transversals;
}

PermSet BSGS::stabilizers(unsigned i) const
{
  return schreier_structures[i]->labels();
}

std::pair<Perm, unsigned> BSGS::strip(Perm const &perm, unsigned offs) const
{
  Perm result(perm);

  for (unsigned i = offs; i < base.size(); ++i) {
    unsigned beta = result[base[i]];
    if (!schreier_structures[i]->contains(beta))
      return std::make_pair(result, i + 1u);

    result *= ~schreier_structures[i]->transversal(beta);
  }

  return std::make_pair(result, base.size() + 1u);
}

bool BSGS::strips_completely(Perm const &perm) const
{
  auto strip_result(strip(perm));

  return strip_result.first.id() && strip_result.second == base.size() + 1u;
}

void BSGS::extend_base(unsigned bp)
{
  base.push_back(bp);

  std::shared_ptr<SchreierStructure> st;

  switch (_transversals) {
    case BSGS::TRANSVERSALS_EXPLICIT:
      st = std::make_shared<ExplicitTransversals>(degree());
      break;
    case BSGS::TRANSVERSALS_SCHREIER_TREES:
    case BSGS::TRANSVERSALS_AUTO:
      st = std::make_shared<SchreierTree>(degree());
      break;
    case BSGS::TRANSVERSALS_SHALLOW_SCHREIER_TREES:
      throw std::logic_error("TODO");
  }

  st->create_root(bp);

  schreier_structures.push_back(st);
}

void BSGS::update_schreier_structure(unsigned i, PermSet const &generators)
{
  unsigned bp = base[i];
  auto st = schreier_structures[i];

  st->create_root(bp);
  st->create_labels(generators);

  std::vector<unsigned> stack {bp};
  std::set<unsigned> done {bp};

  while (!stack.empty()) {
    unsigned beta = stack.back();
    stack.pop_back();

    for (auto i = 0u; i < generators.size(); ++i) {
      Perm gen = generators[i];
      unsigned beta_prime = gen[beta];

      if (done.find(beta_prime) == done.end()) {
        done.insert(beta_prime);
        stack.push_back(beta_prime);

        st->create_edge(beta_prime, beta, i);
      }
    }
  }
}

// TODO: add option to keep original generators
void BSGS::remove_generators()
{
  Dbg(Dbg::DBG) << "Removing redundant strong generators from BSGS:";
  Dbg(Dbg::DBG) << *this;

  Dbg(Dbg::TRACE) << "Stabilizers are:";
#ifndef NDEBUG
  for (auto i = 0u; i < base.size(); ++i)
    Dbg(Dbg::TRACE) << "S(" << i + 1u << ") = " << stabilizers(i);
#endif

  assert(base.size() > 0u);

  unsigned n = strong_generators[0].degree();
  unsigned k = base.size();

  // TODO
  std::unordered_set<Perm> strong_generator_set(
    strong_generators.begin(), strong_generators.end());

  auto difference = [&](std::unordered_set<Perm> const &lhs,
                        std::unordered_set<Perm> const &rhs) {

    int delta = static_cast<int>(lhs.size() - rhs.size());
    std::unordered_set<Perm> res;

    for (auto const &perm : lhs) {
      if (rhs.find(perm) == rhs.end() &&
          strong_generator_set.find(perm) != strong_generator_set.end()) {

        res.insert(perm);
        if (--delta == 0)
          break;
      }
    }

    return res;
  };

  auto produces_orbit = [&](unsigned root,
                            std::unordered_set<Perm> const &generators,
                            std::vector<unsigned> const &orbit) {

    std::vector<int> in_orbit_ref(n + 1u, 0), in_orbit(n + 1u, 0);

    for (unsigned x : orbit)
      in_orbit_ref[x] = 1;

    if (!in_orbit_ref[root])
      return false;

    in_orbit[root] = 1;

    std::vector<unsigned> queue {root};
    auto target = orbit.size() - 1u;

    while (!queue.empty()) {
      unsigned x = queue.back();
      queue.pop_back();

      for (auto const &gen : generators) {
        unsigned y = gen[x];
        if (!in_orbit_ref[y])
          return false;

        if (in_orbit[y] == 0) {
          in_orbit[y] = 1;
          queue.push_back(y);

          if (--target == 0u) {
            return true;
          }
        }
      }
    }

    return false;
  };

  std::unordered_set<Perm> stabilizer_set;
  std::unordered_set<Perm> stabilizer_intersection;

  for (int i = static_cast<int>(k - 1u); i >= 0; --i) {
    auto tmp(stabilizers(i));
    std::unordered_set<Perm> stabilizer_set_next(tmp.begin(), tmp.end());

    stabilizer_intersection = difference(stabilizer_set_next, stabilizer_set);

    stabilizer_set = stabilizer_set_next;

    Dbg(Dbg::TRACE) << "=== Considering S(" << i + 1u << ")/S(" << i + 2u << ")"
                    << " = " << stabilizer_intersection;

    if (stabilizer_intersection.size() < 2u)
      continue;

    auto it(stabilizer_intersection.begin());
    while (it != stabilizer_intersection.end()) {
      Dbg(Dbg::TRACE) << "Considering " << *it;

#ifndef NDEBUG
      std::unordered_set<Perm> reduced_stabilizers(stabilizer_set);

      auto it_(reduced_stabilizers.begin());
      while (it_ != reduced_stabilizers.end()) {
        if (*it_ == *it) {
          reduced_stabilizers.erase(it_++);
          break;
        } else {
          ++it_;
        }
      }
#endif
      auto orbit_gens(stabilizer_set);
      orbit_gens.erase(*it);

      bool remove_stab = false;
      if (produces_orbit(base[i], orbit_gens, orbit(i))) {
#ifndef NDEBUG
        Dbg(Dbg::TRACE) << base[i] << "^" << reduced_stabilizers
                        << " = " << orbit(i);
#endif

        Dbg(Dbg::TRACE) << "=> Removing strong generator " << *it;

        remove_stab = true;
      }
#ifndef NDEBUG
      else {
        Dbg(Dbg::TRACE) << base[i] << "^" << reduced_stabilizers
                        << " =/= " << orbit(i);
      }
#endif

      if (remove_stab) {
        strong_generator_set.erase(*it);
        stabilizer_set.erase(*it);
#ifndef NDEBUG
        reduced_stabilizers.erase(*it);
#endif
        stabilizer_intersection.erase(it++);
      } else {
        ++it;
      }
    }
  }

  strong_generators = PermSet(strong_generator_set.begin(),
                              strong_generator_set.end());

  Dbg(Dbg::TRACE) << "Remaining strong generators: " << strong_generators;

  std::vector<int> stabilizes(strong_generators.size(), 0);
  PermSet stabilizers;

  for (int i = static_cast<int>(base.size()) - 1; i >= 0; --i) {
    for (auto j = 0u; j < strong_generators.size(); ++j) {
      if (stabilizes[j])
        continue;

      bool stab = true;
      for (int k = 0; k < i; ++k) {
        if (strong_generators[j][base[k]] != base[k]) {
          stab = false;
          break;
        }
      }

      if (stab) {
        stabilizes[j] = 1;
        stabilizers.insert(strong_generators[j]);
      }
    }

    Dbg(Dbg::TRACE) << "=> S(" << i + 1u << ") = " << stabilizers;

    update_schreier_structure(i, stabilizers);
  }
}

std::ostream& operator<<(std::ostream& stream, BSGS const &bsgs)
{
  stream << "BASE: [";

  for (auto i = 0u; i < bsgs.base.size(); ++i) {
    stream << bsgs.base[i];
    if (i < bsgs.base.size() - 1u)
      stream << ", ";
  }

  stream << "]; SGS: [";

  for (auto i = 0u; i < bsgs.strong_generators.size(); ++i) {
    stream << bsgs.strong_generators[i];
    if (i < bsgs.strong_generators.size() - 1u)
      stream << ", ";
  }

  stream << ']';
  return stream;
}

} // namespace cgtl
