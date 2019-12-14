#ifndef _GUARD_PROFILE_RUN_H
#define _GUARD_PROFILE_RUN_H

#include <cstdlib>
#include <stdexcept>
#include <string>

#include "profile_timer.h"
#include "profile_util.h"

std::string run_gap(std::string const &script, double *t);

template<typename FUNC>
void run_cpp(FUNC &&f, double *t)
{
  // run group creation in child process
  pid_t maybe_child;
  switch ((maybe_child = timer_start())) {
  case -1:
    throw std::runtime_error("failed to fork child process");
  case 0:
    {
      f();

      _Exit(EXIT_SUCCESS);
    }
    break;
  }

  *t = timer_stop(maybe_child);
}

#endif // _GUARD_PROFILE_RUN_H
