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

  struct S {
    S() {};
    S(std::vector<PartialPerm> const &generators);

    unsigned dom;
    std::vector<unsigned> alpha;
    std::vector<PartialPerm> generators;
    std::vector<std::vector<unsigned>> action_component;
    SchreierTree schreier_tree;
    OrbitGraph orbit_graph;
    std::vector<unsigned> scc;
  };

  class RClass {
  public:
    RClass(PartialPerm const &x, EEMP::S const &s);

    bool is_member(PartialPerm const &s);

  private:
    PartialPerm _x;
    EEMP::S _s;
    std::vector<std::vector<unsigned>> _action_component;
    EEMP::SchreierTree _schreier_tree;
    EEMP::OrbitGraph _orbit_graph;
    PermGroup _schreier_generators;
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

  static std::vector<PartialPerm> inverse_semigroup_r_class_representatives(
    S const &s);

  static bool is_inverse_semigroup_member(PartialPerm const &y,
    EEMP::S const &s, PermGroup const &sx);
};

std::ostream& operator<<(
  std::ostream& stream, EEMP::SchreierTree const &schreier_tree);

std::ostream& operator<<(
  std::ostream& stream, EEMP::OrbitGraph const &orbit_graph);

} // namespace cgtl

#endif // _GUARD_EEMP_H
