#ifndef GUARD_UTIL_H
#define GUARD_UTIL_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

namespace mpsym
{

namespace internal
{

namespace util
{

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

template<typename IT>
std::size_t container_hash(IT first, IT last) {
  // see: https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
  // boost container hash API is not stable across a sufficient range of versions
  std::size_t seed = std::distance(first, last);;
  for (auto it = first; it != last; ++it)
    seed ^= *it + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  return seed;
}

template<typename T>
struct ContainerHash
{
  std::size_t operator()(T const &c) const
  { return util::container_hash(c.begin(), c.end()); }
};

inline std::mt19937 random_engine()
{ return std::mt19937{std::random_device{}()}; }

} // namespace util

} // namespace internal

} // namespace mpsym

#endif // GUARD_UTIL_H
