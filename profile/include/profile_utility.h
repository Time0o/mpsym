#ifndef _GUARD_PROFILE_UTILITY_H
#define _GUARD_PROFILE_UTILITY_H

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include <sys/types.h>

template<typename ...ARGS>
void print(std::ostream &os, char const *prefix, char const *endl, ARGS &&...args)
{
  os << prefix;

  using expander = int[];
  (void)expander{0, (void(os << " " << std::forward<ARGS>(args)), 0)...};

  os << endl;
}

template<typename... ARGS>
void info(ARGS &&...args)
{ print(std::cout, "INFO:", "\n", std::forward<ARGS>(args)...); }

template<typename... ARGS>
void progress(ARGS &&...args)
{ print(std::cout, "INFO:", "\r", std::forward<ARGS>(args)...); }

inline void progress_done()
{ std::cout << "\n"; }

template<typename... ARGS>
void result(ARGS &&...args)
{
  std::cout << std::scientific;
  std::cout.precision(3);

  print(std::cout, "RESULT:", "\n", std::forward<ARGS>(args)...);
}

template<typename... ARGS>
void warning(ARGS &&...args)
{ print(std::cerr, "WARNING:", "\n", std::forward<ARGS>(args)...); }

template<typename... ARGS>
void error(ARGS &&...args)
{ print(std::cerr, "ERROR:", "\n", std::forward<ARGS>(args)...); }

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
