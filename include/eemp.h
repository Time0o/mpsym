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

  static std::pair<unsigned, std::vector<unsigned>>
  strongly_connected_components(OrbitGraph const &orbit_graph);

  static SchreierTree scc_spanning_tree(
    unsigned i, OrbitGraph const &orbit_graph,
    std::vector<unsigned> const &scc);

  static PartialPerm schreier_trace(
    unsigned x, SchreierTree const &schreier_tree,
    std::vector<PartialPerm> const &generators, unsigned dom_max,
    unsigned target = 0u);

  static PermGroup schreier_generators(unsigned i,
    std::vector<PartialPerm> const &generators, unsigned dom_max,
    std::vector<std::vector<unsigned>> const &action_component,
    SchreierTree const &schreier_tree, OrbitGraph const &orbit_graph,
    std::vector<unsigned> const &sccs);

  static std::vector<PartialPerm> r_class_representatives(
    SchreierTree const &schreier_tree,
    std::vector<PartialPerm> const &generators);
};

std::ostream& operator<<(
  std::ostream& stream, EEMP::SchreierTree const &schreier_tree);

std::ostream& operator<<(
  std::ostream& stream, EEMP::OrbitGraph const &orbit_graph);

} // namespace cgtl

#endif // _GUARD_EEMP_H
