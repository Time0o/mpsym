#ifndef GUARD_EEMP_H
#define GUARD_EEMP_H

#include <ostream>
#include <utility>
#include <vector>

#include "partial_perm.hpp"
#include "partial_perm_set.hpp"

namespace mpsym
{

namespace internal
{

class PartialPerm;
class PermGroup;

namespace eemp
{

using Node = std::vector<unsigned>;

using Component = std::vector<Node>;

struct SchreierTree
{ std::vector<std::pair<unsigned, unsigned>> data; };

struct OrbitGraph
{ std::vector<std::vector<unsigned>> data; };

struct Sccs
{
  std::vector<unsigned> data;
  std::vector<std::vector<unsigned>> data_expanded() const;
};

Component action_component(Node const &alpha,
                           PartialPermSet const &generators,
                           SchreierTree &schreier_tree,
                           OrbitGraph &orbit_graph);

Sccs strongly_connected_components(OrbitGraph const &orbit_graph);

// TODO: parameters???
SchreierTree spanning_tree(OrbitGraph const &orbit_graph,
                           Sccs const &sccs,
                           unsigned scc);

PartialPerm schreier_trace(PartialPermSet const &generators,
                           SchreierTree const &schreier_tree,
                           unsigned from,
                           unsigned to);

//PermGroup schreier_generators(
//  unsigned i, std::vector<PartialPerm> const &generators, unsigned dom_max,
//  std::vector<std::vector<unsigned>> const &action_component,
//  SchreierTree const &schreier_tree, OrbitGraph const &orbit_graph,
//  std::vector<unsigned> const &sccs);
//
//std::vector<PartialPerm> r_class_representatives(
//  SchreierTree const &schreier_tree,
//  std::vector<PartialPerm> const &generators);

std::ostream &operator<<(std::ostream &os, SchreierTree const &schreier_tree);

std::ostream &operator<<(std::ostream &os, OrbitGraph const &orbit_graph);

std::ostream &operator<<(std::ostream &os, Sccs const &sccs);

} // namespace eemp

} // namespace internal

} // namespace mpsym

#endif // GUARD_EEMP_H
