#ifndef _GUARD_PROFILE_RUN_H
#define _GUARD_PROFILE_RUN_H

#include <cstdlib>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "profile_timer.h"
#include "profile_util.h"

std::string run_gap(std::string const &script,
                    bool hide_output,
                    bool hide_errors,
                    double *t);

// TODO: use invoke_result instead
template<typename FUNC>
typename std::result_of<FUNC()>::type run_cpp(FUNC &&f, double *t)
{
  // true if this function does not return anything
  constexpr bool returns_void =
    std::is_same<typename std::result_of<FUNC()>::type, void>::value;

  timer_start();

  if (timer_realtime_enabled()) {
    // run group creation in this process
    if constexpr (returns_void) {
      f();

      if (t)
        *t = timer_stop();
    } else {
      auto ret = f();

      if (t)
        *t = timer_stop();

      return ret;
    }
  } else {
    // run group creation in child process
    if constexpr (!returns_void)
      throw std::logic_error("non void functions require realtime timer");

    pid_t child;
    switch ((child = fork())) {
    case -1:
      throw std::runtime_error("failed to fork child process");
    case 0:
      {
        f();

        _Exit(EXIT_SUCCESS);
      }
      break;
    }

    if (t)
      *t = timer_stop(child);
  }
}

#endif // _GUARD_PROFILE_RUN_H
