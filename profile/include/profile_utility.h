#ifndef _GUARD_PROFILE_UTILITY_H
#define _GUARD_PROFILE_UTILITY_H

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include <sys/types.h>

template<typename Arg, typename... Args>
void error(Arg &&arg, Args&&... args)
{
  // TODO: progname

  std::cerr << "error: ";

  std::cerr << std::forward<Arg>(arg);

  using expander = int[];
  (void)expander{0, (void(std::cerr << ' ' << std::forward<Args>(args)), 0)...};

  std::cerr << '\n';
}

template<typename T>
T stox(std::string const &str)
{
  long long i;

  bool success = true;

  try {
    std::size_t idx;
    i = std::stoll(str, &idx);

    success = idx == str.size() &&
              i >= static_cast<long long>(std::numeric_limits<T>::min()) &&
              i <= static_cast<long long>(std::numeric_limits<T>::max());

  } catch (...) {
    success = false;
  }

  if (!success)
    throw std::invalid_argument("stox failed");

  return static_cast<T>(i);
}

#endif // _GUARD_PROFILE_UTILITY_H
