#ifndef _GUARD_EEMP_H
#define _GUARD_EEMP_H

#include <ostream>
#include <utility>
#include <vector>

#include "partial_perm.h"
#include "perm_group.h"

namespace cgtl
{

class EEMP
{
public:
  struct SchreierTree {
    std::vector<std::pair<unsigned, unsigned>> data;
  };

  struct OrbitGraph {
    std::vector<std::vector<unsigned>> data;
  };

  static std::vector<std::vector<unsigned>> action_component(
    std::vector<unsigned> const &alpha,
    std::vector<PartialPerm> const &generators, unsigned dom_max,
    SchreierTree &schreier_tree, OrbitGraph &orbit_graph);

  static std::vector<unsigned> strongly_connected_components(
    OrbitGraph const &orbit_graph);

  static PartialPerm schreier_trace(
    unsigned x, SchreierTree const &schreier_tree,
    std::vector<PartialPerm> const &generators, unsigned dom_max);

  static PermGroup schreier_generators(PartialPerm const &x,
    std::vector<PartialPerm> const &generators, unsigned dom_max,
    std::vector<std::vector<unsigned>> const &action_component,
    SchreierTree const &schreier_tree, OrbitGraph const &orbit_graph);

  static std::vector<PartialPerm> r_class_elements(PartialPerm const &x,
    std::vector<PartialPerm> const &generators, unsigned dom_max);

  static std::vector<PartialPerm> r_classes_in_d_class(PartialPerm const &x,
    std::vector<PartialPerm> const &generators, unsigned dom_max);
};

std::ostream& operator<<(
  std::ostream& stream, EEMP::SchreierTree const &schreier_tree);

std::ostream& operator<<(
  std::ostream& stream, EEMP::OrbitGraph const &orbit_graph);

} // namespace cgtl

#endif // _GUARD_EEMP_H
