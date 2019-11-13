#ifndef _GUARD_PROFILE_RUN_H
#define _GUARD_PROFILE_RUN_H

#include <cstdlib>
#include <string>

#include "profile_timer.h"
#include "profile_utility.h"

bool run_gap(std::string const &script, double *t);

template<typename FUNC>
bool run_cpp(FUNC &&f, double *t)
{
  // run group creation in child process
  pid_t maybe_child;
  switch ((maybe_child = timer_start())) {
  case -1:
    error("failed to fork child process");
    return false;
    break;
  case 0:
    {
      f();

      _Exit(EXIT_SUCCESS);
    }
    break;
  }

  *t = timer_stop(maybe_child);

  return true;
}

#endif // _GUARD_PROFILE_RUN_H
