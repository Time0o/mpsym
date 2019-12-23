#include <chrono>
#include <cstdlib>
#include <ostream>
#include <sstream>
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
  : _enabled(false),
    _timer("profile")
  {}

  void enable()
  { _enabled = true; }

  bool enabled()
  { return _enabled; }

  void start()
  { _timer.start(); }

  void stop(double *t)
  {
    _timer.stop();
    *t = _timer.total();
  }

private:
  bool _enabled;
  timer::Timer _timer;

} realtime_timer;

bool stop_child(pid_t child)
{
  int status;
  waitpid(child, &status, 0);

  return  WIFEXITED(status) && WEXITSTATUS(status) == EXIT_SUCCESS;
}

} // namespace

void timer_realtime_enable()
{ realtime_timer.enable(); }

bool timer_realtime_enabled()
{ return realtime_timer.enabled(); }

void timer_start()
{
  if (realtime_timer.enabled())
    realtime_timer.start();
  else
    child_timer.start();
}

double timer_stop(pid_t child)
{
  double t;
  if (realtime_timer.enabled()) {
    realtime_timer.stop(&t);
  } else {
    if (child == 0)
      throw std::logic_error("no child pid given");

    if (!stop_child(child))
      throw std::runtime_error("the forked child process terminated prematurely");

    child_timer.stop(&t);
  }

  return t;
}

void debug_timer_dump(char const *timer)
{
  std::ostream *os = TIMER_GET_OUT();

  std::stringstream ss;
  TIMER_SET_OUT(&ss);

  TIMER_DUMP(timer);
  debug(ss.str());

  TIMER_SET_OUT(os);
}
