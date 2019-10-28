#ifndef _GUARD_PERM_H
#define _GUARD_PERM_H

#include <cstddef>
#include <functional>
#include <ostream>
#include <vector>

/**
 * @file perm.h
 * @brief Defines `Perm`.
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

class Perm;
class PermWord;

} // namespace cgtl

namespace std
{

template<>
struct hash<cgtl::Perm>
{
  std::size_t operator()(cgtl::Perm const &perm) const;
};

template<>
struct hash<cgtl::PermWord>
{
  std::size_t operator()(cgtl::PermWord const &perm_word) const;
};

} // namespace std

namespace cgtl
{

/** A permutation representation.
 *
 * The `Perm` class represents a permutation on the set \f$X = \{1, \dots,
 * n\}\f$ of positive integers (with an arbitrary fixed \f$n\f$), i.e. a
 * bijective mapping from \f$X\f$ to itself. It realizes common operations like
 * permutation application and chaining via convenient operator overloads.
 * Only permutations on continuous ranges of positive integers starting from
 * \f$1\f$ are supported, i.e. a permutation represented by a `Perm` object
 * always has a domain of the form \f$\{1, 2, 3, \dots, n\}\f$ with \f$ n \in
 * \mathbb{N}_+\f$.
 */
class Perm
{
friend std::size_t std::hash<Perm>::operator()(Perm const &perm) const;

public:
  /** Construct an *identity permutation*.
   *
   * An identity permutation maps every element in its domain to itself.
   *
   * \param degree
   *	 determines the domain of the resulting permutation (which is \f$\{1,
   *	 \dots, degree\}\f$). For `degree == 0u` the permutation' domain is the
   *	 empty set (i.e. using `operator[]` with any index is invalid in this
   *	 case and multiplication with any other permutation will yield another
   *	 permutation over the empty set)
   */
  Perm(unsigned degree = 1);

  /** Construct a permutation from an explicit mapping description.
   *
   * \param perm
   *     a vector of element images, if the image at the \f$i\f$th index has
   *     the value \f$j\f$, the resulting permutation will map \f$i + 1\f$ to
   *     \f$j\f$; all element images must thus be greater than `0u` and
   *     different from each other in order to describe a valid permutation,
   *     otherwise this constructor's behaviour is undefined
   */
  Perm(std::vector<unsigned> const &perm);

  /** Construct a permutation from a permutation description given as product of
   *  cycles.
   *
   * This constructor operates in accordance with the commonly encountered
   * notation of permutations as products of disjoint cycles in mathematical
   * literature (although here the cycles need not be disjoint, in this case
   * the cycles are "chained together left to right" in accordance with
   * `operator*=`). Note that this is mostly a convenience constructor which
   * will almost always be slower than an equivalent call to
   * Perm(std::vector<unsigned> const &perm).
   *
   * \param degree
   *	 determines the domain of the resulting permutation (\f$\{1, \dots,
   *	 degree\}\f$), which can not be unambiguously determined from the given
   *	 `cycles`
   *
   * \param cycles
   *     a cycle description of the resulting permutation, e.g. `{{1, 2}, {3,
   *     4}}` to describe the permutation \f$(1 2)(3 4)\f$ (assuming `degree ==
   *     4u`); if any cycle element is greater than `degree` this constructor's
   *     behaviour is undefined
   */
  Perm(unsigned degree, std::vector<std::vector<unsigned>> const &cycles);

  /** Apply a permutation.
   *
   * The result of invoking `operator[x]` on a permutation \f$p\f$ is the
   * positive integer \f$y\f$ such that \f$y = p(x)\f$. If \f$x \notin dom(p)\f$
   * this function's behaviour is undefined.
   *
   * \param x the positive integer on which to apply this permutation
   *
   * \return the positive integer `y` to which `x` is mapped by this permutation
   */
  unsigned const& operator[](unsigned const x) const;

  /** Construct the inverse of a permutation.
   *
   * A permutation \f$p\f$'s inverse \f$p^{-1}\f$ is the permutation which maps
   * every element \f$y \in im(p)\f$ to \f$x \in dom(p)\f$ if \f$y = p(x)\f$.
   *
   * \return a newly constructed `Perm` object inverse to this one
   */
  Perm operator~() const;

  /** Check whether two permutations are equal.
   *
   * Two permutations are equal by definition if they map all elements in their
   * domains (which must be equal) to the same elements.
   *
   * \param rhs
   *     the permutation which is to be compared with this one, if `rhs`'s
   *     degree() is not equal to this permutation's degree(), this operator's
   *     behaviour is undefined
   *
   * \return `true` if `rhs` is equal to this permutation according to the
   *         definition above
   */
  bool operator==(Perm const &rhs) const;

  /** Check whether two permutations are not equal.
   *
   * This function's result is always equal to `!(*this == rhs)`.
   *
   * \param rhs
   *     the permutation which is to be compared with this one, if `rhs`'s
   *     degree() is not equal to this permutation's degree(), this operator's
   *     behaviour is undefined
   *
   * \return `true` if `rhs` is not equal (see `operator==`) to this permutation
   */
  bool operator!=(Perm const &rhs) const;

  /** "Chain" another permutation to an existing permutation.
   *
   * After invoking `p1 *= p2` for `Perm` objects `p1` and `p2`, `p1` is the
   * permutation which maps every element \f$x \in dom(p_1)\f$ to
   * \f$p2(p1(x))\f$.
   *
   * \param rhs the permutation which is to be "chained" to this one
   *
   * \return a reference to this permutation
   */
  Perm& operator*=(Perm const &rhs);

  /** Obtain a permutation's *degree*.
   *
   * Assuming a permutation's domain is the set \f$\{1, \dots, n\}\f$, its
   * degree is the positive integer \f$n\f$
   *
   * \return this permutation's degree
   */
  unsigned degree() const { return _n; }

  /** Check whether a permutation is an *identity permutation*.
   *
   * An identity permutation maps every element in its domain to itself.
   *
   * \return `true` if this permutation is an identity, else `false`
   */
  bool id() const;

  /** Obtain a permutation's vector representation.
   *
   * \return vector representation of this permutation.
   */
  std::vector<unsigned> vect() const { return _perm; }

  /** Obtain a permutation's cycle representation.
   *
   * \return a cycle representation of this permutation.
   */
  std::vector<std::vector<unsigned>> cycles() const;

  /** *Extend* a permutation's domain.
   *
   * This function constructs a new permutation \f$p'\f$ from this permutation,
   * \f$p\f$, with the domain \f$\{1, \dots, n'\}\f$ where \f$n\f$ and \f$n'\f$
   * are the degrees of \f$p\f$ and \f$p'\f$ with \f$n' \geq n\f$ and the
   * action of \f$p'\f$ is defined by \f$p'(i) = p(i)\f$ for \f$1 \leq i \leq
   * n\f$ and \f$p'(i) = i\f$ for \f$n + 1 \leq i \leq n'\f$.
   *
   * \param degree
   *     degree of the resulting permutation as described above; if degree is
   *     smaller than this permutation's degree, this function's behaviour is
   *     undefined
   *
   * \return a copy of this permutation extended to the domain \f$\{1, \dots,
   *         n'\}\f$
   */
  Perm extended(unsigned degree) const;

  /** *Normalize* a permutation's domain.
   *
   * This function constructs a new permutation \f$p'\f$ from this permutation,
   * \f$p\f$, with the domain \f$\{1, \dots, high - low + 1\}\f$ whose action
   * is defined by \f$p'(i) = j \iff p(i + low) = j + low\f$ for \f$1 \leq i
   * \leq high - low + 1\f$.
   *
   * \param low
   *     the low end of the domain range which is to be "shifted downwards" as
   *     described above
   *
   * \param high
   *     the high end of the domain range which is to be "shifted downwards" as
   *     described above
   *
   * \return a *normalized* copy of this permutation as described above
   */
  Perm normalized(unsigned low, unsigned high) const;

  /** *Shift* a permutation's domain.
   *
   * This function constructs a new permutation \f$p'\f$ equivalent to this
   * permutation, \f$p\f$, but with the domain \f$\{1 + shift, \dots, shift + n -
   * 1\}\f$ (where \f$n\f$ is this permutation's degree).
   *
   * \param shift domain shift as described above
   *
   * \return a *shifted* copy of this permutation as described above
   */
  Perm shifted(unsigned shift) const;

  /** Construct a *restricted* version of a permutation.
   *
   * "Restricting" a permutation to a new domain \f$d_{new}\f$ conceptually has
   * the same effect as removing all cycles from the permutation's disjoint
   * cycle representation which contain at least one element \f$x \notin
   * d_{new}\f$.
   *
   * \param domain
   *     the domain to which to restrict the resulting permutation, must not
   *     contain any element outside this permutation's degree or this
   *     function's behaviour is undefined
   *
   * \return a restricted version of this permutation according to the
   *         definition above
   */
  Perm restricted(std::vector<unsigned> const &domain) const;

private:
  unsigned _n;
  std::vector<unsigned> _perm;
};

class PermWord
{
friend std::size_t std::hash<PermWord>::operator()(PermWord const &permWord) const;

public:
  PermWord(Perm const &perm);
  PermWord(unsigned degree = 1) : PermWord(Perm(degree)) {};
  PermWord(std::vector<unsigned> const &perm) : PermWord(Perm(perm)) {};
  PermWord(unsigned n, std::vector<std::vector<unsigned>> const &cycles)
    : PermWord(Perm(n, cycles)) {};

  unsigned operator[](unsigned const i) const;
  PermWord operator~() const;
  bool operator==(PermWord const &rhs) const;
  bool operator!=(PermWord const &rhs) const;
  PermWord& operator*=(PermWord const &rhs);

  unsigned degree() const { return _n; }
  bool id() const;
  Perm perm() const;

private:
  unsigned _n;
  std::vector<Perm> _perms;
  std::vector<Perm> _invperms;
};

std::ostream& operator<<(std::ostream& stream, Perm const &perm);
Perm operator*(Perm const &lhs, Perm const &rhs);

std::ostream& operator<<(std::ostream& stream, PermWord const &pw);
PermWord operator*(PermWord const &lhs, PermWord const &rhs);

inline std::size_t hash_value(Perm const &perm) { return std::hash<Perm>()(perm); }

} // namespace cgtl

#endif // _GUARD_PERM_H
