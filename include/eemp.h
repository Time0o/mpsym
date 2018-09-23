#ifndef _GUARD_EEMP_H
#define _GUARD_EEMP_H

#include <ostream>
#include <utility>
#include <vector>

#include "partial_perm.h"
#include "perm_group.h"

/**
 * @file eemp.h
 * @brief Definitions related to inverse semigroups of partial permutations.
 *
 * This file defines several auxiliary functions used behind the scenes of
 * PartialPermInverseSemigroup. All functions and their implementations are
 * based on \cite east16.
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

namespace eemp
{

/** Schreier tree data structure.
 *
 * Describes an OrbitGraph spanning tree. Similar to, but not to be confused
 * with cgtl::SchreierTree. Should be treated as opaque by
 * functions other than those defined in eemp.h.
 */
struct SchreierTree {
  std::vector<std::pair<unsigned, unsigned>> data;
};

/** Orbit graph data structure.
 *
 * Describes an *orbit graph* according to the description of
 * action_component(). Should be treated as opaque by functions other than
 * those defined in eemp.h.
 */
struct OrbitGraph {
  std::vector<std::vector<unsigned>> data;
};

/** Compute the *component of the action* of a set of partial permutations (the
 * *generators*) on a set of elements \f$\alpha\f$ as well an associated *orbit
 * graph* and spanning  *Schreier tree*.
 *
 * The action component results from the repeated application of the generators
 * to \f$\alpha\f$ (*application* refers to the operation performed by the
 * PartialPerm::image() function) and all resulting sets of elements until a
 * fixpoint is reached. The resulting orbit graph is the unique directed graph
 * connecting nodes corresponding to the resulting sets with edges labelled
 * with the generators that produce the destination node's set when applied to
 * the source node's set. The Schreier tree is a spanning tree through the
 * orbit graph rooted at the node corresponding to \f$\alpha\f$ (in most cases
 * there are several possible such spanning trees, this function makes no
 * guarantees about which one is returned). Orbit graph and Schreier tree can
 * be treated as opaque by code which uses the functions defined in this header.
 *
 * \param alpha
 *     a set of elements in the form of a vector, the first element in the
 *     resulting action component and root of the resulting orbit graph and
 *     Schreier tree
 *
 * \param generators a set of partial permutations in the form of a vector
 *
 * \param dom_max
 *     the maximum over the values `PartialPerm::dom_max()` of all generators
 *     (passed as an argument to avoid frequent recalculation)
 *
 * \param[out] schreier_tree the resulting Schreier tree
 *
 * \param[out] orbit_graph the resulting orbit graph
 *
 * \return the resulting action component in form of a vector of vectors;
 *         Schreier tree and orbit graph are internally represented using
 *         indices into this vector and are thus only meaningful in combination
 *         with it
 */
std::vector<std::vector<unsigned>> action_component(
  std::vector<unsigned> const &alpha,
  std::vector<PartialPerm> const &generators, unsigned dom_max,
  SchreierTree &schreier_tree, OrbitGraph &orbit_graph);

/** Find the *strongly connected components* of an orbit graph.
 *
 * The strongly connected components (short s.c.c.'s) are a partition of the
 * orbit graph such that all nodes in a strongly connected components can be
 * directly or indirectly reached from every other node in the strongly
 * connected component.
 *
 * \param orbit_graph orbit graph for which the s.c.c.'s are to be determined
 *
 * \return a pair in which the first element is the number of resulting s.c.c's
 *         and the second element is a vector containing as many elements as
 *         there are nodes in the orbit graph (and thus elements in the
 *         corresponding action component) in which two elements have the same
 *         value if and only if the nodes corresponding to their indices in the
 *         vector lie in the same s.c.c.
 *
 * Example:
 *
 * Assume that an action component \f$\{\alpha_1, \alpha_2, \alpha_3,
 * \alpha_4\}\f$ is partitioned according to the s.c.c.'s of an associated
 * orbit graph as follows: \f$\{\{\alpha_1\}, \{\alpha_2, \alpha_3\},
 * \{\alpha_4\}\}\f$. This function might then return `{3u, {0u, 1u, 1u, 2u}}`
 * (the concrete values assigned to each s.c.c. in the resulting vector are not
 * guaranteed).
 */
std::pair<unsigned, std::vector<unsigned>> strongly_connected_components(
  OrbitGraph const &orbit_graph);

/** Find a spanning tree for a strongly connected component inside an orbit
 *  graph.
 *
 * \param i
 *     the resulting spanning tree will be calculated for the s.c.c. with index
 *     `i` (i.e. all the elements at indices corresponding to elements in the
 *     s.c.c. given by the `scc` argument have value `i`)
 *
 * \param orbit_graph the orbit graph
 *
 * \param scc
 *     the orbit graphs s.c.c.'s, determined via strongly_connected_components()
 *
 * \return a spanning Schreier tree for the given s.c.c. rooted at the node
 *         inside the s.c.c. which corresponds to the element with the smallest
 *         index in the associated action component
 */
SchreierTree scc_spanning_tree(
  unsigned i, OrbitGraph const &orbit_graph,
  std::vector<unsigned> const &scc);

/** Trace a Schreier tree
 *
 * Trace a Schreier tree for some orbit graph from a node back to its root,
 * i.e. determine the partial permutation which results from chaining the
 * partial permuation edge labels on the way through the Schreier tree from the
 * root to the node together.
 *
 * \param x
 *     the action component index of the node from which the backtracing
 *     operation through the Schreier tree should be performed.
 *
 * \param schreier_tree the Schreier tree
 *
 * \param generators
 *     the generators which were used in the generation of the Schreier (using
 *     action_component()) and form the Schreier tree's edge labels
 *
 * \param dom_max
 *     the maximum over the values `PartialPerm::dom_max()` of all generators
 *     (passed as an argument to avoid frequent recalculation)
 *
 * \param target
 *     the index at which the backtracing operation should terminate, i.e. the
 *     root node's index. Ordinarily this should be `0u` except when this
 *     function is used to perform backtracing through s.c.c. spanning trees
 *     obtain via scc_spanning_tree()
 */
PartialPerm schreier_trace(
  unsigned x, SchreierTree const &schreier_tree,
  std::vector<PartialPerm> const &generators, unsigned dom_max,
  unsigned target = 0u);

/** Compute the *Schreier generators* for one strongly connected component of an
 *  orbit graph.
 *
 * The Schreier generators form a permutation group, for details on their
 * theoretical significance see \cite east16.
 *
 * \param i
 *     the resulting Schreier generators will be calculated for the s.c.c with
 *     index `i` (i.e. all the elements at indices corresponding to elements in
 *     the s.c.c. given by the `scc` argument have value `i`)
 *
 * \param generators
 *     the generators used to generate the orbit graph and the corresponding
 *     action component
 *
 * \param dom_max
 *     the maximum over the values `PartialPerm::dom_max()` of all generators
 *     (passed as an argument to avoid frequent recalculation)
 *
 * \param action_component the action component associated with the orbit graph
 *
 * \param schreier_tree
 *     a spanning Schreier tree through the s.c.c. rooted at the s.c.c.'s element
 *     with the smallest action component index
 *
 * \param orbit_graph the orbit graph
 *
 * \param sccs
 *     the orbit graphs s.c.c.'s, determined via strongly_connected_components()
 */
PermGroup schreier_generators(
  unsigned i, std::vector<PartialPerm> const &generators, unsigned dom_max,
  std::vector<std::vector<unsigned>> const &action_component,
  SchreierTree const &schreier_tree, OrbitGraph const &orbit_graph,
  std::vector<unsigned> const &sccs);

/** Calculate the \f$\mathscr{R}\f$-class respresentatives of an inverse partial
 *  permutation semigroup from a corresponding Schreier tree.
 *
 * See \cite east16 for details on the theory begin \f$\mathscr{R}\f$-classes.
 * In the special case of an inverse semigroup of partial permutations \f$G\f$
 * the \f$\mathscr{R}\f$-class representatives can be determined by first
 * obtaining a Schreier tree for on orbit graph constructed using a generating
 * set for the group and rooted at \f$\cup_{g \in G} dom(g)\f$ and then
 * calculating the partial permutations that result if for every path from the
 * Schreier tree's root to one its leaves the corresponding partial permutation
 * edge labels are chained together.
 *
 * \param schreier_tree
 *     the Schreier tree rooted at \f$\cup_{g \in G} dom(g)\f$ which is traced
 *     in order to determine the \f$\mathscr{R}\f$-class representatives
 *
 * \param generators
 *     the partial permutations which constitute the Schreier tree's edge labels
 *
 * \return a vector containing the \f$\mathscr{R}\f$-class representatives in no
 *         particular order
 */
std::vector<PartialPerm> r_class_representatives(
  SchreierTree const &schreier_tree,
  std::vector<PartialPerm> const &generators);

/** Print a Schreier tree according to the notation in \cite east16.
 *
 * \param stream a stream object
 *
 * \param schreier_tree the Schreier tree
 *
 * \return `stream`
 */
std::ostream& operator<<(
  std::ostream& stream, SchreierTree const &schreier_tree);

/** Print an orbit graph according to the notation in \cite east16.
 *
 * \param stream a stream object
 *
 * \param orbit_graph the orbit graph
 *
 * \return a reference to `stream`
 */
std::ostream& operator<<(
  std::ostream& stream, OrbitGraph const &orbit_graph);

} // namespace eemp

} // namespace cgtl

#endif // _GUARD_EEMP_H
