#ifndef _GUARD_EEMP_H
#define _GUARD_EEMP_H

#include <ostream>
#include <utility>
#include <vector>

#include "partial_perm.h"

namespace cgtl
{

class EEMP
{
public:
  struct SchreierTree {
    std::vector<std::pair<unsigned, unsigned>> data;
    unsigned dom_max;
  };

  struct OrbitGraph {
    std::vector<std::vector<unsigned>> data;
  };

  static std::vector<std::vector<unsigned>> action_component(
    std::vector<unsigned> const &alpha,
    std::vector<PartialPerm> const &generators,
    SchreierTree &schreier_tree, OrbitGraph &orbit_graph);

  static std::vector<unsigned> strongly_connected_components(
    OrbitGraph const &orbit_graph);

  static PartialPerm schreier_trace(
    unsigned x, SchreierTree const &schreier_tree,
    std::vector<PartialPerm> const &generators);

  static std::vector<PartialPerm> schreier_generators(
    PartialPerm const &x, std::vector<PartialPerm> const &generators);
};

std::ostream& operator<<(
  std::ostream& stream, EEMP::SchreierTree const &schreier_tree);

std::ostream& operator<<(
  std::ostream& stream, EEMP::OrbitGraph const &orbit_graph);

} // namespace cgtl

#endif // _GUARD_EEMP_H
