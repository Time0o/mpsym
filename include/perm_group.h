#ifndef _GUARD_PERM_GROUP_H
#define _GUARD_PERM_GROUP_H

#include <map>
#include <vector>

#include "bsgs.h"
#include "perm.h"
#include "schreier_sims.h"

/**
 * @file perm_group.h
 * @brief Defines `PermGroup`.
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

/** A permutation group representation.
 *
 * This class provides a useful abstraction encapsulating several complex
 * algorithms and data structures used to efficiently represent a permutation
 * group defined by a set of generating permutations without the need to store
 * elements explicitely for very large groups.
 */
class PermGroup
{
public:
  enum ConstructionMethod {
    SCHREIER_SIMS,
    SCHREIER_SIMS_RANDOM,
    AUTO
  };

  class const_iterator
  {
  public:
    const_iterator() : _end(true) {};
    const_iterator(PermGroup const &pg);

    const_iterator operator++();
    const_iterator operator++(int) { next_state(); return *this; }
    Perm const & operator*() const { return _current_result; }
    Perm const * operator->() const { return &_current_result; }
    bool operator==(const_iterator const &rhs) const;
    bool operator!=(const_iterator const &rhs) const { return !((*this) == rhs); }

  private:
    void next_state();
    void update_result();

    std::vector<unsigned> _state;
    bool _trivial;
    bool _end;

    std::vector<std::vector<Perm>> _transversals;
    std::vector<Perm> _current_factors;
    Perm _current_result;
  };

  /** Construct the permutation group \f$\{(1)\}\f$.
   */
  PermGroup() : PermGroup(1, {}) {}

  /** Construct a permutation group.
   *
   * Constructs a permutation group representation from a given set of
   * generating permutations. The generators and group elements might not be
   * stored explicitly in the resulting object. Instead some variation of the
   * *Schreier-Sims* algorithm (see \cite holt05, chapter 4.4.2) might be used
   * to compute a *base* and *strong generating* set for the group which
   * describes the group completely and can be used to, among others, test
   * element membership and iterate through all group elements efficiently.
   *
   * \param degree the permutation group's *degree*, must be the same as the
                   degree of all given generators (which in turn implies that
                   they map from the set \f$\{1, \dots, degree\}\f$ to itself)
                   otherwise this constructor's behaviour is undefined
   * \param generators a generating set for the permutation group
   *
   * \param schreier_sims_method
   *     determines which variant of the Schreier-Sims algorithms is used to
   *     construct the internal group representation, might influence this
   *     constructor's runtime
   */
  PermGroup(unsigned degree, std::vector<Perm> const &generators,
    ConstructionMethod = AUTO);

  /** Check two permutation groups for equality.
   *
   * Two permutation groups are equal exactly when they contain the same
   * elements. Although it it possible to use the PermGroup class to represent
   * permutation groups containing a very large number of elements, this
   * operation is guaranteed to be performed in \f$O(|bsgs.sgs()|)\f$ time. If
   * `(*this).degree() != rhs.degree()`, this function's behaviour is undefined.
   *
   * \return `true`, if `*this` and `rhs` describe permutation groups containing
   *         the same elements, else `false`
   */
  bool operator==(PermGroup const &rhs) const;

  /** Check two permutation groups for inequality.
   *
   * The result is always equal to `!(*this == rhs)`.
   *
   * \return `true`, if `*this` and `rhs` describe permutation groups not
   *         containing the same elements, else `false`
   */
  bool operator!=(PermGroup const &rhs) const;

  /** Construct a symmetric permutation group
   *
   * \param degree
   *     degree \f$n\f$ of the resulting group, for `degree == 0u` this
   *     functions behaviour is undefined
   *
   * \return the symmetric group \f$S_n\f$
   */
  static PermGroup symmetric(unsigned degree);

  /** Construct a cyclic permutation group
   *
   * \param degree
   *     degree \f$n\f$ of the resulting group, for `degree == 0u` this
   *     functions behaviour is undefined
   *
   * \return the cyclic group \f$C_n\f$
   */
  static PermGroup cyclic(unsigned degree);

  /** Construct an alternating permutation group
   *
   * \param degree
   *     degree \f$n\f$ of the resulting group, for `degree < 3u` this functions
   *     behaviour is undefined
   *
   * \return the alternating group \f$A_n\f$
   */
  static PermGroup alternating(unsigned degree);

  /** Obtain a constant iterator iterating over this group's elements.
   *
   * Note that the permutation group elements might not be stored explicitly
   * and could instead be constructed by the iterator "on the fly" which means
   * that storing references or pointers to the partial permutation pointed to
   * by any iterator iterating over a permutation group will result in
   * undefined behaviour. The elements are not guaranteed to be returned in any
   * particular order but every element will be returned exactly once.The
   * elements will also be returned in the same order on subsequent iterations
   * through the permutation group. This last point is not necessarily true for
   * other equal (in the sense of `operator==`) permutation groups.
   *
   * Iterating through a permutation group could thus look like this:
   *
   * ~~~~~{.cpp}
   * PermGroup pg(degree, generators);
   *
   * for (auto const &perm : pg)
   *   std::cout << perm;
   *
   * for (auto const it = pg.begin(); it != pg.end(); ++i)
   *   std::cout << *it;
   * ~~~~~
   *
   * \return a constant iterator pointing to some element in this permutation
   *         group, incrementing it will yield all other group members exactly
   *         once (in no particular order) until end() is reached.
   */
  const_iterator begin() const { return const_iterator(*this); }

  /** Obtain a contant iterator signifying a permutation group's "end".
   *
   * \return a contant iterator signifying a permutation group's "end"
   */
  const_iterator end() const { return const_iterator(); }

  /** Obtain a permutation group's *degree*.
   *
   * A permutation group's degree is the positive integer \f$n\f$ such that
   * all its elements map from the set \f$\{1, \dots, n\}\f$ to itself.
   *
   * \return this permutation group's degree
   */
  unsigned degree() const { return _n; }

  /** Obtain a permutation group's *order*.
   *
   * A permutation group \f$G\f$ 's order is equal to the number of elements it
   * contains, i.e. \f$|G|\f$. Note that every permutation group must at least
   * contain an identity and thus it is always true that `order() > 0u`.
   *
   * \return this permutation group's order
   */
  unsigned order() const { return _order; }

  /** Check whether a permutation group is *trivial*.
   *
   * A permutation group is trivial by definition if it only contains an
   * identity permutation on the set \f$\{1, \dots, n\}\f$ (where \f$n\f$ is the
   * group's degree()). This is equivalent to testing whether `order() == 1u`.
   *
   * \return `true` if this permutation group is trivial, else `false`
   */
  bool is_trivial() const { return _bsgs.strong_generators.empty(); }

  /** Check whether a permutation group is symmetric.
   *
   * \return `true` if the permutation group \f$G \leq Sym(\Omega)\f$
   *          represented by this object is in fact equal to \f$Sym(\Omega)\f$,
   *          else `false`
   */
  bool is_symmetric() const;

  /** Check whether a permutation group is symmetric.
   *
   * \return `true` if the permutation group \f$G \leq S_n\f$ represented by
   *          this object is the alternating group \f$A_n\f$, else `false`
   */
  bool is_alternating() const;

  /** Check whether a permutation group is *transitive*.
   *
   * A permutation group is transitive by definition if for the group orbit of
   * any of its elements \f$x\f$ it holds that: \f$G(x) = \{g \cdot x \in \{1,
   * \dots, n\} : g \in G\} = \{1, \dots, n\}\f$ (where \f$n\f$ is the group's
   * degree()).
   *
   * \return `true` if this permutation group is transitive, else `false`.
   */
  bool is_transitive() const;

  /** Obtain a permutation group's base and strong generating set.
   *
   * This function is only meaningful is the permutation group's elements are
   * not stored explicitely, which can be determined via TODO.
   *
   * \return a BSGS object representing this permutation group's base and strong
   *         generating set
   */
  BSGS bsgs() const { return _bsgs; }

  /** Return a permutation group's orbit partition.
   *
   * The *orbit* of a group element \f$g \in G\f$ where \f$G\f$ acts on the set
   * \f$\Omega\f$ is the set \f$G(x) = \{x^g \mid g \in G\}\f$. This function
   * returns the set \f$\{\{x^g \mid g \in G\} \mid x \in \Omega\}\f$ which
   * contains every possible element orbit and is always a partition of \f$G\f$.
   *
   * \return this permutation group's orbit partition as a vector of vectors
   *         where each element vector contains all elements in one of the
   *         possible element orbits (sorted in ascending order) and the element
   *         vectors themselves are sorted by their the smallest element (in
   *         ascending order)
   */
  std::vector<std::vector<unsigned>> orbits() const;

  /** Check whether a permutation group contains a given permutation.
   *
   * Note that the group's elements may not be stored explicitly so while
   * efficient, this operation is not trivial. if `perm.degree() !=
   * (*this).degree()` this function's behaviour is undefined.   *
   * \return `true` if this permutation group contains the permutation `perm`,
   *         else `false`
   */
  bool contains_element(Perm const &perm) const;

  /** Construct a random group element.
   *
   * This function makes no guarantees about the returned elements distribution
   * of cryptographic security. Repeated calls to this function may not be very
   * efficient, use TODO instead.
   *
   * \return some element \f$x\f$ of this permutation group
   */
  Perm random_element() const;

  /** Find a *disjoint subgroup decomposition* of a permutation group
   *
   * This function's implementation is based on \cite donaldson09.
   *
   * A disjoint subgroup decomposition of a permutation group \f$G\f$ is a set
   * \f$\{H_1, \dots, H_k\}\f$ where \f$H_1, \dots, H_k\f$ are subgroups of
   * \f$G\f$ such that \f$G = \{\alpha_1 \cdot \alpha_2 \dots \alpha_k :
   * \alpha_i \in H_i, (1 \leq i \leq k)\}\f$ and \f$\forall (i, j) \in \{1,
   * \dots, k\}^2 : i \neq j \Rightarrow moved(H_i) \cap moved(H_j) =
   * \emptyset\f$.  Here, \f$moved(H_i)\f$ is the set \f$\{x : (x \in \{1,
   * \dots, n\} \land \exists \alpha \in H_i : \alpha(x) \neq x)\}\f$.
   *
   * \param complete
   *    if this is `true` the function will always return the *finest* possible
   *    disjoint subgroup decomposition (i.e. one in which no subgroup can be
   *    further decomposed into disjoint subgroups), this might be more
   *    computationally expensive
   *
   * \param disjoint_orbit_optimization
   *    if this is `true`, the optimization described in \cite donaldson09,
   *    chapter 5.2.1 is applied, this is only meaningful if `complete = true`,
   *    otherwise this argument is ignored.
   *
   * \return a disjoint subgroup decomposition of this permutation group given
   *         as a vector of subgroups as described above, in no particular
   *         order, if no disjoint subgroup partition could be found, the vector
   *         contains the group itself as its only element.
   */
  std::vector<PermGroup> disjoint_decomposition(
    bool complete = true, bool disjoint_orbit_optimization = false) const;

  /** Find a *wreath product decomposition* of a permutation group
   *
   * This function's implementation is based on \cite donaldson09.
   *
   * A wreath product decomposition of a permutation group \f$G \leq S_n\f$ is
   * (using the notation in \cite donaldson09) a triple \f$(H, K,
   * \mathcal{X})\f$. Here\f$H \leq S_m\f$, \f$K \leq S_d\f$ and \f$n = m \cdot
   * d\f$, \f$\mathcal{X}\f$ is a set \f$\{X_1, \dots, X_d\}\f$ of sets
   * \f$X_i\f$ (with \f$|X_i| = m\f$) which partition the set \f$X = \{1,
   * \dots, n\}\f$ on which \f$G\f$ operates and \f$G = \{\sigma(\beta)\
   * \sigma_1(\alpha_1) \dots \sigma_d(\alpha_d) : \beta \in K, \alpha_i \in H
   * \ (1 \leq i \leq d)\}\f$. \f$\sigma\f$ is the permutation representation
   * of the action of \f$K\f$ on \f$\mathcal{X}\f$ which permutes the sets
   * \f$X_i\f$ "among themselves" and the \f$\sigma_i\f$ are the permutation
   * representations of the obvious actions of \f$K\f$ on the sets \f$X_i\f$.
   * This is conventionally written as: \f$G = H \wr K\f$.
   *
   * \return a wreath product decomposition of \f$G\f$ as described above in
   *         the form of the vector of permutation groups \f$[\sigma(K),
   *         \sigma_1(H_1), \dots, \sigma_d(H_d)]\f$, if no wreath product
   *         decomposition could be found (either because none exists or the
   *         algorithm is unable to determine a valid decomposition, which is
   *         currently possible due to limitations in the implementation) an
   *         empty vector is returned
   */
  std::vector<PermGroup> wreath_decomposition() const;

private:
  std::vector<PermGroup> disjoint_decomposition_complete(
    bool disjoint_orbit_optimization = true) const;

  std::vector<PermGroup> disjoint_decomposition_incomplete() const;

  unsigned _n;
  unsigned _order;
  BSGS _bsgs;
};

} // namespace cgtl

#endif // _GUARD_PERM_GROUP_H
