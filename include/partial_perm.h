#ifndef _GUARD_PARTIAL_PERM_H
#define _GUARD_PARTIAL_PERM_H

#include <cstddef>
#include <ostream>
#include <vector>

#include "perm.h"

/**
 * @file partial_perm.h
 * @author Timo Nicolai
 *
 * @brief Defines `PartialPerm`.
 */

namespace cgtl
{

class PartialPerm;

} // namespace cgtl

namespace std
{

/** A partial permutation hash functor.
 */
template<>
struct hash<cgtl::PartialPerm>
{
  std::size_t operator()(cgtl::PartialPerm const &pperm) const;
};

} // namespace std

namespace cgtl
{

class PartialPerm
{
friend std::size_t std::hash<PartialPerm>::operator()(
  PartialPerm const &perm) const;

public:
  PartialPerm(unsigned degree = 0);
  PartialPerm(std::vector<unsigned> const &dom,
              std::vector<unsigned> const &im);
  PartialPerm(std::vector<unsigned> const &pperm);

  static PartialPerm id(std::vector<unsigned> const &dom);
  static PartialPerm from_perm(Perm const &perm);

  unsigned operator[](unsigned const i) const;
  PartialPerm operator~() const;
  bool operator==(PartialPerm const &rhs) const;
  bool operator!=(PartialPerm const &rhs) const;
  PartialPerm& operator*=(PartialPerm const &rhs);

  /** Obtain a partial permutation's domain.
   *
   * \return an ordered vector containing all elements in this partial
   *         permutation's domain in ascending order
   * \sa dom_min, dom_max
   */
  std::vector<unsigned> dom() const { return _dom; }

  /** Obtain the smallest element in a partial permutation's domain.
   *
   * \return the smallest element in this partial permutation's domain or `0u`
             if this partial permutation is empty
   * \sa dom, dom_max
   */
  unsigned dom_min() const { return _dom_min; }

  /** Obtain the largest element in a partial permutation's domain.
   *
   * \return the largest element in this partial permutation' domain or `0u` if
             this partial permutation is empty
   * \sa dom, dom_min
   */
  unsigned dom_max() const { return _dom_max; }

  /** Obtain a partial permutation's image.
   *
   * \return an ordered vector containing all elements in this partial
   *         permutation's image in ascending order
   * \sa im_min, im_max
   */
  std::vector<unsigned> im() const { return _im; }

  /** Obtain the smallest element in a partial permutation's image.
   *
   * \return the smallest element in this partial permutation's image or `0u` if
             this partial permutation is empty
   * \sa im, im_max
   */
  unsigned im_min() const { return _im_min; }

  /** Obtain the largest element in a partial permutation's image.
   *
   * \return the largest element in this partial permutation's image or `0u` if
             this partial permutation is empty
   * \sa im, im_min
   */
  unsigned im_max() const { return _im_max; }

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
  bool empty() const { return _pperm.empty(); }

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
  bool id() const { return _id; }

  /** Contruct a *restricted* version of a partial permutation.
   *
   * *Restriction* in this context means that the resulting partial
   * permutation's domain will be the intersection of this partial
   * permutation's domain and the `domain` parameter (which need not be
   * ordered). The mapping from domain to image elements within this
   * intersection is preserved.
   *
   * Example:
   *
   * ~~~~~{.cpp}
   * PartialPerm pp({1, 2, 3}, {4, 7, 6});
   * PartialPerm pp_restricted(pp.restricted({1, 3})); // result is [1 4][3 6]
   * ~~~~~
   *
   * \param domain vector of elements of domain to which the resulting partial
   *               permutation is to be restricted
   * \return a newly constructed partial permutation equal to this partial
   *         permutation restricted to `domain`
   */
  PartialPerm restricted(std::vector<unsigned> const &domain) const;

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

  /** Compute the image of a vector under this partial permutation.
   *
   * The partial permutation is, if applicable, applied to every element in
   * `alpha` (in order), i.e. for every element \f$x \in \alpha\f$, if \f$x \in
   * dom(*this)\f$, \f$(*this)(x)\f$ is appended to the result.
   */
  std::vector<unsigned> image(std::vector<unsigned> const &alpha) const;

private:
  void update_limits();

  std::vector<unsigned> _pperm;
  std::vector<unsigned> _dom, _im;
  unsigned _dom_min, _dom_max;
  unsigned _im_min, _im_max;
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
std::ostream& operator<<(std::ostream& stream, PartialPerm const &pperm);

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

} // namespace cgtl

#endif // _GUARD_PARTIAL_PERM_H
