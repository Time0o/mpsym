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
  struct SchreierTree {
    unsigned dom_max;
    std::vector<std::pair<unsigned, unsigned>> data;
  };

  struct OrbitGraph {
    unsigned dom_max;
    std::vector<std::vector<unsigned>> data;
  };

  static std::vector<std::vector<unsigned>> action_components(
    std::vector<unsigned> const &alpha,
    std::vector<PartialPerm> const &generators,
    SchreierTree &schreier_tree, OrbitGraph &orbit_graph);

  static PartialPerm schreier_trace(
    SchreierTree const &schreier_tree, unsigned i,
    std::vector<PartialPerm> const &generators);

  static std::vector<std::vector<unsigned>> strongly_connected_components(
    OrbitGraph const &orbit_graph);
};

} // namespace cgtl

#endif // _GUARD_EEMP_H
