#ifndef _GUARD_UTIL_H
#define _GUARD_UTIL_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

/**
 * @file util.h
 * @brief Defines several utility functions used throughout this project.
 *
 * @author Timo Nicolai
 */

namespace util
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

    assert(std::numeric_limits<T>::max() / base >= res);

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
  while (x > 1u) {
    assert(std::numeric_limits<T>::max() / x >= res);

    res *= x--;
  }

  return res;
}


template<typename T, typename U = double>
void mean_stddev(std::vector<T> const &vals, U *mean, U *stddev) {
  T zero = static_cast<T>(0);
  T size = static_cast<T>(vals.size());

  *mean = std::accumulate(vals.begin(), vals.end(), zero) / size;

  std::vector<U> d(vals.size());
  std::transform(vals.begin(), vals.end(), d.begin(),
                 [mean](U val) { return val - *mean; });

  *stddev = std::sqrt(
    std::inner_product(d.begin(), d.end(), d.begin(), zero) / size);
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
