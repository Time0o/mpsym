#ifndef _GUARD_UTIL_H
#define _GUARD_UTIL_H

#include <algorithm>
#include <cstddef>
#include <vector>

/**
 * @file util.h
 * @brief Defines several utility functions used throughout this project.
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

/** Calculate integer powers.
 *
 * `T` must be a signed or unsigned integer type, otherwise this function's
 * behaviour is undefined.
 *
 * \param base base, for negative bases this function's behaviour is undefined
 *
 * \param exp
 *     exponent, for negative exponents or when `base^exp` is not representable
 *     by `T` this function's behaviour is undefined
 *
 * \return base^exp
 */
template<typename T>
T pow(T base, T exp)
{
  T res = 1;

  for (;;) {
    if (exp & 1)
      res *= base;

    if (!(exp >>= 1))
      break;

    base *= base;
  }

  return res;
}

/** Calculate the factorial of a non-negative integer.
 *
 * `T` must be a signed or unsigned integer type, otherwise this function's
 * behaviour is undefined.
 *
 * \param x
 *     the non-negative input number \f$x\f$, if \f$x < 0\f$ or \f$x!\f$ is not
 *     representable by `T`, $this function's behaviour is undefined
 *
 * \return \f$x!\f$
 */
template<typename T>
T factorial(T x)
{
  T res = 1u;
  while (x > 1u)
    res *= x--;

  return res;
}

/** Expand a compact set partition representation into an explicit which allows
 *  iteration over the partitions.
 *
 * T must be an unsigned integer type, `I` must by a type implementing
 * `operator==` (usually also some integer type).
 *
 * \param partition
 *     a partition in form of a vector in which set elements correspond to
 *     indices and elements belonging to the same partition are marked by equal
 *     vector elements at their repective indices
 *
 * \return a vector of vectors, where each subvector contains all subelements in
 *         a single partition (in ascending order) and in which the subvectors
 *         are ordered according to their minimum elements (in ascending order)
 */
template<typename T, typename I>
std::vector<std::vector<T>> expand_partition(std::vector<I> partition)
{
  auto mm = std::minmax_element(partition.begin(), partition.end());
  auto first = *mm.first;
  auto last = *mm.second;

  std::vector<std::vector<T>> res(last - first + 1);
  for (T x = 0u; x < partition.size(); ++x)
    res[partition[x] - first].push_back(x);

  std::sort(res.begin(), res.end(),
           [](std::vector<T> a, std::vector<T> b){ return a[0] < b[0]; });

  return res;
}

/** Compute a hash value for a vector of (hashable) elements.
 *
 * T must be a hashable type, the result depends both on the values contained in
 * `vec` and their order.
 *
 * \param vec a vector of elements
 *
 * \return a sufficiently good hash over `vec`
 */
template<typename T>
std::size_t vector_hash(std::vector<T> const &vec) {
  // see: https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
  // boost container hash API is not stable across a sufficient range of versions
  std::size_t seed = vec.size();
  for (auto const &x : vec)
    seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  return seed;
}

} // namespace cgtl

#endif // _GUARD_UTIL_H
