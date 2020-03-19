#include <chrono>
#include <cstdlib>
#include <stdexcept>

#include <sys/times.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "profile_timer.h"
#include "profile_util.h"
#include "timer.h"

namespace
{

class ChildTimer
{
public:
  void start() const
  {}

  void stop(double *t)
  {
    struct tms tms;
    times(&tms);

    *t = static_cast<double>(tms.tms_cutime) / sysconf(_SC_CLK_TCK) - _acc;
    _acc += *t;
  }

private:
  double _acc = 0.0;

} child_timer;

class RealtimeTimer
{
public:
  RealtimeTimer()
  : _timer("profile")
  {}

  void start()
  { _timer.start(); }

  void stop(double *t)
  {
    _timer.stop();
    *t = _timer.last();
  }

private:
  timer::Timer _timer;

} realtime_timer;

} // anonymous namespace

namespace profile
{

void timer_start()
{
#ifdef PROFILE_CPU_TIMER
  child_timer.start();
#else
  realtime_timer.start();
#endif // PROFILE_CPU_TIMER
}

double timer_stop(pid_t child)
{
  double t;

#ifdef PROFILE_CPU_TIMER
  if (child == 0)
    throw std::logic_error("no child pid given");

  int status;
  waitpid(child, &status, 0);

  bool child_stopped =  IFEXITED(status) && WEXITSTATUS(status) == EXIT_SUCCESS;

  if (!child_stopped)
    throw std::runtime_error("the forked child process terminated prematurely");

  child_timer.stop(&t);
#else
  (void)child;

  realtime_timer.stop(&t);
#endif // PROFILE_CPU_TIMER

  return t;
}

} // namespace profile
