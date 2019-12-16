#include <algorithm>
#include <cassert>
#include <unordered_set>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "orbits.h"
#include "perm.h"
#include "perm_set.h"

namespace cgtl
{

void BSGS::reduce_gens()
{
  DBG(DEBUG) << "Removing redundant strong generators from BSGS:";
  DBG(DEBUG) << *this;

  DBG(TRACE) << "Stabilizers are:";
#ifndef NDEBUG
  for (auto i = 0u; i < base_size(); ++i)
    DBG(TRACE) << "S(" << i + 1u << ") = " << stabilizers(i);
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

    DBG(TRACE) << "=== Considering S(" << i + 1u << ")/S(" << i + 2u << ")"
               << " = " << stabilizer_intersection;

    if (stabilizer_intersection.size() < 2u)
      continue;

    auto it(stabilizer_intersection.begin());
    while (it != stabilizer_intersection.end()) {
      DBG(TRACE) << "Considering " << *it;

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
      PermSet orbit_gens(stabilizer_set.begin(), stabilizer_set.end());
      orbit_gens.erase(*it);

      bool remove_stab = false;
      if (orbit_check(base_point(i), orbit_gens, orbit(i))) {
#ifndef NDEBUG
        DBG(TRACE) << base_point(i) << "^" << reduced_stabilizers
                   << " = " << orbit(i);
#endif

        DBG(TRACE) << "=> Removing strong generator " << *it;

        remove_stab = true;
      }
#ifndef NDEBUG
      else {
        DBG(TRACE) << base_point(i) << "^" << reduced_stabilizers
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

  DBG(DEBUG) << "Reduced BSGS:";
  DBG(DEBUG) << *this;
}

std::unordered_set<Perm> BSGS::reduce_gens_set_difference(
  std::unordered_set<Perm> const &lhs,
  std::unordered_set<Perm> const &rhs,
  std::unordered_set<Perm> const &base) const
{
  std::unordered_set<Perm> res;

  for (auto const &perm : lhs) {
    if (rhs.find(perm) == rhs.end() && base.find(perm) != base.end())
      res.insert(perm);
  }

  return res;
}


} // namespace cgtl
