#ifndef GUARD_PARTIAL_PERM_H
#define GUARD_PARTIAL_PERM_H

#include <cstddef>
#include <ostream>
#include <set>
#include <vector>

#include "perm.h"

/**
 * @file partial_perm.h
 * @author Timo Nicolai
 *
 * @brief Defines `PartialPerm`.
 */

namespace mpsym
{

namespace internal
{

class PartialPerm;

} // namespace internal

} // namespace mpsym

namespace std
{

/** A partial permutation hash functor.
 */
template<>
struct hash<mpsym::internal::PartialPerm>
{
  std::size_t operator()(mpsym::internal::PartialPerm const &pperm) const;
};

} // namespace std

namespace mpsym
{

namespace internal
{

class PartialPerm
{
friend std::size_t std::hash<PartialPerm>::operator()(
  PartialPerm const &perm) const;

public:
  PartialPerm(unsigned degree = 0u);
  PartialPerm(std::vector<unsigned> const &dom,
              std::vector<unsigned> const &im);
  PartialPerm(std::vector<unsigned> const &pperm);

  unsigned operator[](unsigned const i) const;
  PartialPerm operator~() const;
  bool operator==(PartialPerm const &rhs) const;
  bool operator!=(PartialPerm const &rhs) const;
  PartialPerm& operator*=(PartialPerm const &rhs);

  static PartialPerm from_perm(Perm const &perm);

  /** Construct a permutation from a partial permutation.
   *
   * This operation is only meaningful if the partial permutation only contains
   * cycles (i.e. does not contain chains). Otherwise, the result is undefined.
   * The resulting permutation contains the same cycles as the partial
   * permutation. (provided they contain only elements \f$x \leq degree\f$)
   * and maps all elements \f$x\f$ with \f$1 < x < degree\f$ and \f$x \notin
   * dom(*this)\f$ to themselves.
   *
   * \param degree degree of the resulting permutation
   * \return permutation with degree `degree` corresponding to this partial
   *         permutation
   */
  Perm to_perm(unsigned degree) const;

  /** Obtain a partial permutation's domain.
   *
   * \return an ordered vector containing all elements in this partial
   *         permutation's domain in ascending order
   * \sa dom_min, dom_max
   */
  std::vector<unsigned> dom() const
  { return _dom; }

  /** Obtain the smallest element in a partial permutation's domain.
   *
   * \return the smallest element in this partial permutation's domain or `0u`
             if this partial permutation is empty
   * \sa dom, dom_max
   */
  unsigned dom_min() const
  { return _dom.empty() ? 0u : _dom[0]; }

  /** Obtain the largest element in a partial permutation's domain.
   *
   * \return the largest element in this partial permutation' domain or `0u` if
             this partial permutation is empty
   * \sa dom, dom_min
   */
  unsigned dom_max() const
  { return _dom.empty() ? 0u : _dom.back(); }

  /** Obtain a partial permutation's image.
   *
   * \return an ordered vector containing all elements in this partial
   *         permutation's image in ascending order
   * \sa im_min, im_max
   */
  std::vector<unsigned> im() const
  { return _im; }

  /** Obtain the smallest element in a partial permutation's image.
   *
   * \return the smallest element in this partial permutation's image or `0u` if
             this partial permutation is empty
   * \sa im, im_max
   */
  unsigned im_min() const
  { return _im.empty() ? 0u : _im[0]; }

  /** Obtain the largest element in a partial permutation's image.
   *
   * \return the largest element in this partial permutation's image or `0u` if
             this partial permutation is empty
   * \sa im, im_min
   */
  unsigned im_max() const
  { return _im.empty() ? 0u : _im.back(); }

  /** Check whether a partial permutation is *empty*.
   *
   * A partial permutation is empty exactly when it is the unique partial
   * permutation which maps from the empty set to itself. In this case dom()
   * and im() will return empty vectors and dom_min(), dom_max(), im_min() and
   * im_max() will all return `0u`.
   *
   * An empty partial permutation can be constructed in any of the following
   * ways:
   *
   * ~~~~~{.cpp}
   * PartialPerm pp_empty1();
   * PartialPerm pp_empty2(0);
   * PartialPerm pp_empty3({});
   * PartialPerm pp_empty4({}, {});
   * PartialPerm pp_empty5(PartialPerm::id({}));
   * ~~~~~
   *
   * \return `true` if this partial permutation is empty, else `false`
   */
  bool empty() const
  { return _pperm.empty(); }

  /** Check whether a partial permutation is an identity.
   *
   * A partial permutation is an identity exactly when it maps every element
   * in its domain to itself. An empty partial permutation is always an
   * identity. Notice also that unlike in the case of permutations, there is no
   * single unique identity partial permutation.
   *
   * Identity permutations (in this case `(1)(2)(3)`) can be constructed in any
   * of the following ways:
   *
   * ~~~~~{.cpp}
   * PartialPerm pp_id1(3);
   * PartialPerm pp_id4({1, 2, 3});
   * PartialPerm pp_id3({1, 2, 3}, {1, 2, 3});
   * PartialPerm pp_id5(PartialPerm::id({1, 2, 3}));
   * ~~~~~
   *
   * \return `true` if this partial permutation is an identity, else `false`
   */
  bool id() const
  { return _id; }

  // TODO
  template<typename IT>
  PartialPerm restricted(IT first, IT last) const
  {
    if (first == last)
      return PartialPerm();

    std::vector<unsigned> pperm_restricted(dom_max(), 0u);

    for (IT it = first; it != last; ++it) {
      unsigned x = *it;

      if (x < dom_min() || x > dom_max())
        continue;

      unsigned y = (*this)[x];

      if (y != 0u)
        pperm_restricted[x - 1u] = y;
    }

    if (!pperm_restricted.empty()) {
      unsigned dom_max_restricted = dom_max();

      while (pperm_restricted[dom_max_restricted - 1u] == 0u)
        --dom_max_restricted;

      pperm_restricted.resize(dom_max_restricted);
    }

    return PartialPerm(pperm_restricted);
  }

  // TODO
  template<template<typename ...> class T, typename IT>
  T<unsigned> image(IT first, IT last) const // TODO
  {
    std::set<unsigned> res;

    for (IT it = first; it != last; ++it) {
      unsigned x = *it;

      if (x < dom_min() || x > dom_max())
        continue;

      unsigned y = _pperm[x - 1u];

      if (y != 0u)
        res.insert(y);
    }

    return T<unsigned>(res.begin(), res.end());
  }

private:
  std::vector<unsigned> _pperm, _dom, _im;
  bool _id;
};

/** Print a partial permutation.
 *
 * Partial permutations are represented in the chain/cycle notation commonly
 * found in mathematical literature.
 *
 * \param stream a stream object
 *
 * \param pperm a `PartialPerm` object
 *
 * \return a reference to `stream`
 */
std::ostream &operator<<(std::ostream &os, PartialPerm const &pperm);

/** Chain two partial permutations together.
 *
 * \param lhs a `PartialPerm` object representing a partial permuation \f$f\f$
 *
 * \param rhs a `PartialPerm` object representing a partial permuation \f$g\f$
 *
 * \return a `PartialPerm` object representing the partial permuation
 *         \f$g \cdot f\f$ with \f$(g \cdot f)(x) = g(f(x)) \in
 *         im(g|_{im(f)})\f$ for \f$x \in dom(f)\f$
 */
PartialPerm operator*(PartialPerm const &lhs, PartialPerm const &rhs);

} // namespace internal

} // namespace mpsym

#endif // GUARD_PARTIAL_PERM_H
