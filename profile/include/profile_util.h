#ifndef _GUARD_PROFILE_UTIL_H
#define _GUARD_PROFILE_UTIL_H

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <sys/types.h>

namespace profile
{

template<typename ...ARGS>
void print(std::ostream &os, char const *prefix, char const *endl, ARGS &&...args)
{
  os << prefix;

  using expander = int[];
  (void)expander{0, (void(os << " " << std::forward<ARGS>(args)), 0)...};

  os << endl << std::flush;
}

template<typename... ARGS>
void info(ARGS &&...args)
{ print(std::cout, "INFO:", "\n", std::forward<ARGS>(args)...); }

template<typename... ARGS>
void result(ARGS &&...args)
{
  std::cout << std::scientific;
  std::cout.precision(3);

  print(std::cout, "RESULT:", "\n", std::forward<ARGS>(args)...);
}

template<typename... ARGS>
void debug(ARGS &&...args)
{ print(std::cout, "DEBUG:", "\n", std::forward<ARGS>(args)...); }

template<typename... ARGS>
void debug_progress(ARGS &&...args)
{ print(std::cout, "DEBUG:", "\r", std::forward<ARGS>(args)...); }

inline void debug_progress_done()
{ std::cout << "\n"; }

void debug_timer_dump(char const *timer);

template<typename... ARGS>
void warning(ARGS &&...args)
{ print(std::cerr, "WARNING:", "\n", std::forward<ARGS>(args)...); }

template<typename... ARGS>
void error(ARGS &&...args)
{ print(std::cerr, "ERROR:", "\n", std::forward<ARGS>(args)...); }

template<typename T>
T stox(std::string const &str)
{
  T i;
  bool success;

  try {
    std::size_t idx;

    if (std::is_signed<T>::value)
      i = std::stoll(str, &idx);
    else
      i = std::stoull(str, &idx);

    success = idx == str.size();
  } catch (...) {
    success = false;
  }

  if (!success)
    throw std::invalid_argument("stox failed");

  return i;
}

std::vector<std::string> split(std::string const &str, char const *delim = " ");

std::string join(std::vector<std::string> const &strs, char const *delim = ",");

} // namespace profile

#endif // _GUARD_PROFILE_UTIL_H
