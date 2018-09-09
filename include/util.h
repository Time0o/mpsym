#ifndef _GUARD_UTIL_H
#define _GUARD_UTIL_H

#include <algorithm>
#include <cstddef>
#include <vector>

namespace cgtl
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

    base *= base;
  }

  return res;
}

template<typename T>
T factorial(T x)
{
  T res = 1u;
  while (x > 1u)
    res *= x--;

  return res;
}

template<typename T>
std::size_t vector_hash(std::vector<T> const &vec) {
  // see: https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
  // boost container hash API is not stable across a sufficient range of versions
  std::size_t seed = vec.size();
  for (auto const &x : vec)
    seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  return seed;
}

template<typename T>
std::vector<std::vector<T>> expand_partition(std::vector<T> set)
{
  auto mm = std::minmax_element(set.begin(), set.end());
  auto first = *mm.first;
  auto last = *mm.second;

  std::vector<std::vector<T>> res(last - first + 1);
  for (auto x = 0u; x < set.size(); ++x)
    res[set[x] - first].push_back(x);

  std::sort(res.begin(), res.end(),
           [](std::vector<T> a, std::vector<T> b){ return a[0] < b[0]; });

  return res;
}

} // namespace cgtl

#endif // _GUARD_UTIL_H
