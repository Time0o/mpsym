#include <algorithm>
#include <cassert>
#include <unordered_set>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "perm.h"
#include "perm_set.h"

namespace cgtl
{

void BSGS::reduce_gens(bool preserve_original)
{
  if (preserve_original)
    throw std::logic_error("TODO");

  if (base_size() == 0u)
    return; // TODO: remove?

  Dbg(Dbg::DBG) << "Removing redundant strong generators from BSGS:";
  Dbg(Dbg::DBG) << *this;

  Dbg(Dbg::TRACE) << "Stabilizers are:";
#ifndef NDEBUG
  for (auto i = 0u; i < base_size(); ++i)
    Dbg(Dbg::TRACE) << "S(" << i + 1u << ") = " << stabilizers(i);
#endif

  std::unordered_set<Perm> strong_generator_set(_strong_generators.begin(),
                                                _strong_generators.end());

  std::unordered_set<Perm> stabilizer_set;
  std::unordered_set<Perm> stabilizer_intersection;

  for (int i = static_cast<int>(base_size() - 1u); i >= 0; --i) {
    auto stabilizer_vect(stabilizers(i));

    std::unordered_set<Perm> stabilizer_set_next(stabilizer_vect.begin(),
                                                 stabilizer_vect.end());

    stabilizer_intersection = reduce_gens_set_difference(stabilizer_set_next,
                                                         stabilizer_set,
                                                         strong_generator_set);

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
      if (reduce_gens_produces_orbit(base_point(i), orbit_gens, orbit(i))) {
#ifndef NDEBUG
        Dbg(Dbg::TRACE) << base_point(i) << "^" << reduced_stabilizers
                        << " = " << orbit(i);
#endif

        Dbg(Dbg::TRACE) << "=> Removing strong generator " << *it;

        remove_stab = true;
      }
#ifndef NDEBUG
      else {
        Dbg(Dbg::TRACE) << base_point(i) << "^" << reduced_stabilizers
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

  _strong_generators = PermSet(strong_generator_set.begin(),
                               strong_generator_set.end());

  Dbg(Dbg::DBG) << "Reduced BSGS:";
  Dbg(Dbg::DBG) << *this;

  reduce_gens_redetermine_schreier_structures();
}

std::unordered_set<Perm> BSGS::reduce_gens_set_difference(
  std::unordered_set<Perm> const &lhs,
  std::unordered_set<Perm> const &rhs,
  std::unordered_set<Perm> const &base) const
{
  int delta = static_cast<int>(lhs.size() - rhs.size());

  assert(delta > 0);

  std::unordered_set<Perm> res;

  for (auto const &perm : lhs) {
    if (rhs.find(perm) == rhs.end() && base.find(perm) != base.end()) {
      res.insert(perm);
      if (--delta == 0)
        break;
    }
  }

  return res;
}

bool BSGS::reduce_gens_produces_orbit(unsigned root,
                                      std::unordered_set<Perm> const &generators,
                                      std::vector<unsigned> const &orbit) const
{
  std::vector<int> in_orbit_ref(degree() + 1u, 0);
  std::vector<int> in_orbit(degree() + 1u, 0);

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
}

void BSGS::reduce_gens_redetermine_schreier_structures()
{
  Dbg(Dbg::TRACE) << "Updating Schreier structures:";

  std::vector<int> stabilizes(_strong_generators.size(), 0);

  PermSet stabilizers;

  for (int i = static_cast<int>(base_size()) - 1; i >= 0; --i) {
    for (auto j = 0u; j < _strong_generators.size(); ++j) {
      if (stabilizes[j])
        continue;

      bool stab = true;
      for (int k = 0; k < i; ++k) {
        if (_strong_generators[j][base_point(k)] != base_point(k)) {
          stab = false;
          break;
        }
      }

      if (stab) {
        stabilizes[j] = 1;
        stabilizers.insert(_strong_generators[j]);
      }
    }

    Dbg(Dbg::TRACE) << "=> S(" << i + 1u << ") = " << stabilizers;

    update_schreier_structure(i, stabilizers);
  }
}

} // namespace cgtl
