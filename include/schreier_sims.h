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
  enum Variant { SIMPLE, RANDOM };

  static std::vector<unsigned> orbit(unsigned alpha,
    std::vector<Perm> const &generators, SchreierTree &st);

  static std::pair<Perm, unsigned> strip(Perm const &perm,
    std::vector<unsigned> const &base, std::vector<SchreierTree> const &sts);

  static void schreier_sims(std::vector<unsigned> &base,
    std::vector<Perm> &generators, std::vector<SchreierTree> &sts);

  static void schreier_sims_random(std::vector<unsigned> &base,
    std::vector<Perm> &generators, std::vector<SchreierTree> &sts,
    unsigned w = 10);
};

} // namespace cgtl

#endif // _GUARD_SCHREIER_SIMS_H
