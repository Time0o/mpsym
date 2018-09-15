#ifndef _GUARD_SCHREIER_SIMS_H
#define _GUARD_SCHREIER_SIMS_H

#include <map>
#include <utility>
#include <vector>

#include "perm.h"
#include "schreier_sims.h"

/**
 * @file schreier_sims.h
 * @brief Definitions related to the Schreier-Sim algorithm.
 *
 * For details on the Schreier-Sims algorithm consult \cite holt05.
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

namespace schreier_sims
{

/** Data structure for space efficient storage of orbit transversals.
 *
 * The *orbit* of a group element \f$g \in G\f$ where \f$G\f$ acts on the set
 * \f$\Omega\f$ is the set \f$G(x) = \{x^g \mid g \in G\}\f$. The *orbit
 * transversals* corresponding to this orbit are \f$\{u_{\beta} \mid \beta \in
 * G(x)\}\f$ where \f$u_{\beta} \in G\f$ maps \f$x\f$ to \f$\beta\f$. Instead
 * of storing these explicitly (which can might consume a lot of space for
 * large orbits (possible for large sets \f$\Omega\f$), this data structure
 * only stores a spanning tree for the orbit graph connecting every element
 * \f$\beta' \in G(x)\f$ to every other element \f$\beta\f$ via a directed edge
 * labelled with \f$g \in G(x)\f$ if \f$\exists g \in G \colon \beta^g = \beta'\f$.
 * The spanning tree is rooted at \f$x\f$.
 */
struct SchreierTree
{
  SchreierTree(unsigned degree) : _degree(degree) {}

  /** Set a spanning tree's root.
   *
   * \param root
   *     the element \f$x\f$ as per the notation used in this struct's
   *     definition
   */
  void create_root(unsigned root) { _root = root; }

  /** Add an edge to a Schreier tree.
   *
   * \param origin
   *     an element \f$\beta' \in \Omega\f$ as per the notation used in this
   *     struct's definition
   *
   * \param destination
   *     an element \f$\beta \in \Omega\f$ as per the notation used in this
   *     struct's definition
   *
   * \param perm
   *     a perm object describing a permutation \f$g \in G\f$ with \f$\beta^g =
   *     \beta'\f$
   */
  void create_edge(unsigned origin, unsigned destination, Perm const &perm) {
    _edges[origin] = destination;
    _labels[origin] = perm;
  }

  /** Obtain a Schreier tree's root element.
   *
   * \return the Schreier tree's root, i.e. the element \f$x\f$ as per the
   *         notation used in this struct's definition
   */
  unsigned root() const { return _root; }

  /** Obtain a Schreier tree's nodes.
   *
   * \return the Schreier tree's nodes, i.e. the orbit \f$G(x)\f$ as per the
   *         notation used in this struct's definition, in the form of a vector
   *         with no particular guaranteed order
   */
  std::vector<unsigned> nodes() const {
    std::vector<unsigned> result {_root};

    for (auto const &node : _edges)
      result.push_back(node.first);

    return result;
  }

  /** Check whether an element of \f$\Omega\f$ is a node of a Schreier tree.
   *
   * This is equivalent to testing whether the element is contained in the orbit
   * \f$G(x)\f$.
   *
   * \param node node to be testing for membership in this Schreier tree
   *
   * \return `true` if `node` is a node of this Schreier tree, else `false`
   */
  bool contains(unsigned node) const {
    return (node == _root) || (_edges.find(node) != _edges.end());
  }

  /** Obtain an orbit transversal from a Schreier tree.
   *
   * \param origin
   *     an element \f$\beta\f$ of \f$\Omega\f$, this must be a node of this
   *     Schreier tree, otherwise this function's behaviour is undefined
   *
   * \return a permutation \f$g \in G\f$ such that \f$x^g = \beta\f$
   */
  Perm transversal(unsigned origin) const
  {
    Perm result(_degree);

    unsigned current = origin;
    while(current != _root) {
      Perm const &label = _labels.find(current)->second;
      result = label * result;
      current = _edges.find(current)->second;
    }

    return result;
  }

private:
  unsigned _degree;
  unsigned _root = 0;
  std::map<unsigned, unsigned> _edges;
  std::map<unsigned, Perm> _labels;
};

/** *Schreier-Sims* variant constants.
 *
 *  These Constants should be used by functions and classes defined in other
 *  files when differentiating between different variants of the
 *  *Schreier-Sims* algorithm.
 */
enum Variant {
  /** Refers to the simple *Schreier-Sims* algorithm implemented in
   *  schreier_sims()
   */
  SIMPLE,
  /** Refers to the *Random Schreier-Sims* algorithm implemented in
   *  schreier_sims_random()
   */
  RANDOM
};

/** Compute the orbit decomposition of a permutation group given as a set of
 *  generating permutations.
 *
 * The *orbit* of a group element \f$g \in G\f$ where \f$G\f$ acts on the set
 * \f$\Omega\f$ is the set \f$G(x) = \{x^g \mid g \in G\}\f$. This function
 * returns the set \f$\{\{x^g \mid g \in G\} \mid x \in \Omega\}\f$ which
 * contains every possible element orbit and is always a partition of \f$G\f$.
 *
 * \param generators
 *     a vector of `Perm` objects describing a generating set for the
 *     permutation group \f$G\f$
 *
 * \return \f$G\f$'s orbit partition as a vector of vectors where each element
 *         vector contains all elements in one of the possible element orbits
 *         (sorted in ascending order) and the element vectors themselves are
 *         sorted by their the smallest element (in ascending order)
 */
std::vector<std::vector<unsigned>> orbits(std::vector<Perm> const &generators);

/** Compute the orbit of a single element \f$x \in \Omega\f$ for the group \f$G
 *  \leq Sym(\Omega)\f$ described by a set of generating permutations.
 *
 * The *orbit* is the set \f$G(x) = \{x^g \mid g \in G\}\f$. This function also
 * computes the corresponding *orbit transversals*, i.e. the set \f$\{u_{\beta}
 * \mid \beta \in G(x)\}\f$ where \f$u_{\beta} \in G\f$ maps \f$x\f$ to
 * \f$\beta\f$. These orbit transversals are stored in a *Schreier tree*
 * (see schreier_sims() and SchreierTree).
 *
 * \param x
 *     the element \f$x \in \Omega\f$ for which \f$G(x)\f$ should be computed
 *
 * \param generators
 *     a vector of `Perm` objects describing a generating set for the
 *     permutation group \f$G\f$
 *
 * \param[out] st
 *     a pointer to a preallocated Schreier tree which will store the orbit
 *     transversals after the execution of this function, provided
 *     `st != nullptr`
 *
 * \return the orbit \f$G(x)\f$ as described above in the form of an ordered
 *         vector of positive integers \f$\in \Omega\f$
 */
std::vector<unsigned> orbit(
  unsigned x, std::vector<Perm> const &generators, SchreierTree *st = nullptr);

/** Test membership of an arbitrary element in \f$Sym(\Omega)\f$ in a
 *  permutation group \f$G\f$ acting on \f$\Omega\f$.
 *
 * This function computes a pair containing a permutation \f$h\f$ and a
 * positive integer \f$i\f$. If the permutation \f$g\f$ described by `perm` is
 * contained in the group \f$G\f$ with *base* `base` and *base orbits/base
 * orbit transversals* described by the *Schreier trees* `sts` (see also
 * schreier_sims()), \f$h\f$ will be the identity permutation and \f$i\f$ will
 * be \f$k + 1\f$ (where \f$k\f$ is the length of the given base). This implies
 * that \f$g\f$ can be decomposed as \f$u_k u_{k-1} \dots u_i\f$ with \f$u_i
 * \in U^{(i)}\f$ where \f$U^{(i)}\f$ is the \f$i\f$th *base orbit transversal*.
 *
 * If \f$g \notin G\f$, this will not hold. The returned values \f$h\f$ and
 * \f$i\f$ are in this case nevertheless vital to the implementation of the
 * *Schreier-Sims* algorithm in schreier_sims() (This is an implementation
 * detail, refer to \cite holt05 for further information).
 *
 * \param perm
 *     the permutation which should be tested for membership in the group with
 *     base `base` and base orbits/base orbit transversals described by `sts`
 *
 * \param base a group base, refer to schreier_sims()
 *
 * \param sts a vector of Schreier trees, refer to schreier_sims()
 *
 * \return a pair containing a permutation \f$h\f$ and a positive integer
 *         \f$i\f$ as described above
 */
std::pair<Perm, unsigned> strip(
    Perm const &perm, std::vector<unsigned> const &base,
    std::vector<SchreierTree> const &sts);

/** The simple *Schreier-Sims* algorithm
 *
 * The following explanation is based on \cite holt05, consult the book for
 * a more detailed depiction.
 *
 * This algorithm computes a *base and strong generating set* (short *BSGS*)
 * completely representing a specific permutation group.  A base of a
 * permutation group \f$\f$ is a list \f$B = [\beta_1, \dots, \beta_k]\f$ of
 * elements \f$\beta_i \in \Omega,\ 1 \leq i \leq k\f$ (where \f$\Omega\f$ is
 * the set on which \f$G\f$ acts) where the only element in \f$G\f$ that fixes
 * every element in \f$B\f$ is the identity.
 *
 * A strong generating set is a set of partial permutations such that for \f$i
 * = 1, 2, \dots, k\f$ it holds that \f$G^{(i)} = \left<S^{(i)}\right>\f$.
 * Here, \f$G^{(i)}\f$ is the \f$i\f$th *base stabilizer* which maps all
 * elements \f$\beta_1, \dots, \beta_{i-1}\f$ to themselves and \f$S^{(i)} := S
 * \cap G^{(i)}\f$.
 *
 * This function also calculates *base orbits/base orbit transversals* in the form of
 * a vector of Schreier trees (one per base element). The \f$i\f$th base orbit
 * is the set \f$\Delta^{(i)} := \beta_i^{H^{(i)}}\f$ (where \f$H^{(i)} :=
 * \left<S^{(i)}\right>\f$) and the corresponding base orbit transversal is
 * \f$U^{(i)} := \{u_{\beta}^{(i)} \mid \beta \in \Delta^{(i)}\}\f$ (where
 * \f$u_{\beta}^{(i)} \in G\f$ maps \f$\beta_i\f$ to \f$\beta\f$).
 *
 * Note that in general there is no unique BSGS for a given permutation group
 * and this function makes no guarantees as to which possible BSGS is returned.
 *
 * \param[in,out] base
 *     initial base, the resulting base will be an extension of this base, i.e.
 *     in general more eleme elements will be appended to the end of this
 *     parameter; base may be empty (and should be if initially no particular
 *     base prefix is desired)
 *
 * \param[in,out] generators
 *     the generating set for the permutation for which BSGS should be computed.
 *     this function replaces the element in this vector by the determined
 *     strong generators (in no particular order)
 *
 * \param[out] sts
 *     schreier trees which describe the orbits of base elements and the
 *     associated transversals which are needed for several further computations
 */
void schreier_sims(
  std::vector<unsigned> &base, std::vector<Perm> &generators,
  std::vector<SchreierTree> &sts);

/** The *Random Schreier-Sims* algorithm
 *
 * This is a variant of the algorithm implemented by schreier_sims() which uses
 * randomly generated group elements (via PrRandomizer) to obtain a BSGS. The
 * generated BSGS is not guaranteed to be correct BSGS for the permutation group
 * described by the given generators, to achieve certainty, an additional call
 * to TODO is necessary.
 *
 * \param[in,out] base
 *     initial base in general more eleme elements will be appended to the end
 *     of this parameter; base may be empty (and should be if initially no
 *     particular base prefix is desired)
 *
 * \param[in,out] generators
 *     the generating set for the permutation for which BSGS should be computed.
 *     this function replaces the element in this vector by the determined
 *     strong generators (in no particular order)
 *
 * \param[out] sts
 *     schreier trees which describe the orbits of base elements and the
 *     associated transversals which are needed for several further computations
 *
 * \param w
 *     tuning parameter, higher values increase the probability that the
 *     computed BSGS is actually a BSGS for the permutation group generated by
 *     the given generators; more specifically this probability is approximately
 *     \f$1 - 2^{-w}\f$, refer to \cite holt05 for details
 */
void schreier_sims_random(
  std::vector<unsigned> &base, std::vector<Perm> &generators,
  std::vector<SchreierTree> &sts, unsigned w = 10);

} // namespace schreier_sims

} // namespace cgtl

#endif // _GUARD_SCHREIER_SIMS_H
