#ifndef GUARD_NUMERIC_H
#define GUARD_NUMERIC_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace mpsym
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

} // namespace util

} // namespace mpsym

#endif // GUARD_NUMERIC_H
