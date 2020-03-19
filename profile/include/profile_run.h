#ifndef _GUARD_PROFILE_RUN_H
#define _GUARD_PROFILE_RUN_H

#include <cstdlib>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "profile_timer.h"
#include "profile_util.h"

namespace profile
{

std::vector<std::string> run_gap(std::initializer_list<std::string> packages,
                                 std::string const &script,
                                 unsigned num_discarded_runs,
                                 unsigned num_runs,
                                 bool hide_output,
                                 bool hide_errors,
                                 std::vector<double> *ts = nullptr);

// TODO: use invoke_result instead
template<typename FUNC>
typename std::result_of<FUNC()>::type run_cpp(FUNC &&f,
                                              unsigned num_discarded_runs,
                                              unsigned num_runs,
                                              std::vector<double> *ts = nullptr)
{
  // true if this function does not return anything
  constexpr bool returns_void =
    std::is_same<typename std::result_of<FUNC()>::type, void>::value;

  for (unsigned r = 0u; r < num_discarded_runs + num_runs; ++r) {
    bool take_time = ts && r >= num_discarded_runs;

    if (take_time)
      timer_start();

#ifdef PROFILE_CPU_TIMER
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

    if (take_time)
      ts->push_back(timer_stop(child));
  }
#else
    // run group creation in this process
    if constexpr (returns_void) {
      f();

      if (take_time)
        ts->push_back(timer_stop());

    } else {
      auto ret = f();

      if (take_time)
        ts->push_back(timer_stop());

      if (r == num_runs - 1u)
        return ret;
    }
  }

  if constexpr (!returns_void)
    throw std::logic_error("unreachable");
#endif
}

} // namespace profile

#endif // _GUARD_PROFILE_RUN_H
