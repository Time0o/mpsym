#ifndef _GUARD_EEMP_H
#define _GUARD_EEMP_H

#include <utility>
#include <vector>

#include "partial_perm.h"

namespace cgtl
{

class EEMP
{
public:
  static std::vector<std::vector<unsigned>> action_components(
    std::vector<unsigned> const &alpha, std::vector<PartialPerm> const &generators,
    std::vector<std::pair<unsigned, unsigned>> &schreier_tree,
    std::vector<std::vector<unsigned>> &orbit_graph);
};

} // namespace cgtl

#endif // _GUARD_EEMP_H
