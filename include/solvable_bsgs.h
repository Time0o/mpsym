#ifndef _GUARD_SOLVABLE_BSGS_H
#define _GUARD_SOLVABLE_BSGS_H

#include <vector>

#include "perm.h"
#include "schreier_sims.h"

namespace cgtl
{

std::vector<Perm> normalizing_generator(
  Perm const &gen, std::vector<unsigned> &base, std::vector<Perm> &generators,
  std::vector<schreier_sims::SchreierTree> &sts);

} // namespace cgtl

#endif // _GUARD_SOLVABLE_BSGS_H
