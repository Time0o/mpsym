#ifndef _GUARD_PROFILE_UTIL_H
#define _GUARD_PROFILE_UTIL_H

#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <sys/types.h>

#include "timer.hpp"

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

inline void debug_timer_dump(char const *timer)
{
  if (!TIMER_EXISTS(timer)) {
    debug("TIMER (" + std::string(timer) + "): never invoked");
    return;
  }

  std::ostream *os = TIMER_GET_OUT();

  std::stringstream ss;
  TIMER_SET_OUT(&ss);

  TIMER_DUMP(timer);

  auto str(ss.str());
  debug(str.substr(0u, str.size() - 1u));

  TIMER_SET_OUT(os);
}

template<typename... ARGS>
void warning(ARGS &&...args)
{ print(std::cerr, "WARNING:", "\n", std::forward<ARGS>(args)...); }

template<typename... ARGS>
void error(ARGS &&...args)
{ print(std::cerr, "ERROR:", "\n", std::forward<ARGS>(args)...); }

} // namespace profile

#endif // _GUARD_PROFILE_UTIL_H
