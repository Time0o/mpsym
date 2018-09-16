#ifndef _GUARD_BSGS_H
#define _GUARD_BSGS_H

#include <vector>

#include "perm.h"
#include "schreier_sims.h"

namespace cgtl
{

struct BSGS
{
  std::vector<unsigned> base;
  std::vector<Perm> strong_generators;
  std::vector<schreier_sims::SchreierTree> schreier_trees;

  bool contains(Perm const &perm) {
    auto strip_result(schreier_sims::strip(perm, base, schreier_trees));
    return strip_result.first.id() && strip_result.second == base.size() + 1u;
  }

  std::vector<unsigned> orbit(unsigned i) const {
    return schreier_trees[i].nodes();
  }

  Perm transversal(unsigned i, unsigned o) const {
    return schreier_trees[i].transversal(o);
  }

  std::vector<Perm> transversals(unsigned i) const {
    std::vector<Perm> transversals;
    for (unsigned o : orbit(i))
      transversals.push_back(schreier_trees[i].transversal(o));

    return transversals;
  }

  std::vector<Perm> stabilizers(unsigned i) const {
    return schreier_trees[i].labels();
  }
};

} // namespace cgtl

#endif // _GUARD_BSGS_H
