#ifndef _GUARD_UTIL_H
#define _GUARD_UTIL_H

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

} // namespace cgtl

#endif // _GUARD_UTIL_H
