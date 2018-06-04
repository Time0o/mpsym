#ifndef _GUARD_SCHREIER_SIMS_H
#define _GUARD_SCHREIER_SIMS_H

#include <utility>
#include <vector>

#include "perm.h"
#include "schreier_tree.h"

namespace cgtl
{

class SchreierSims
{
public:
  static std::vector<unsigned> orbit(unsigned alpha,
    std::vector<Perm> const &generators, SchreierTree &st);

  static std::pair<Perm, unsigned> strip(Perm const &perm,
    std::vector<unsigned> const &base, std::vector<SchreierTree> const &sts);

  static void schreier_sims(std::vector<unsigned> &base,
    std::vector<Perm> &generators, std::vector<SchreierTree> &sts);
};

} // namespace cgtl

#endif // _GUARD_SCHREIER_SIMS_H
