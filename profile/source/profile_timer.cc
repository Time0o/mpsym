#include <chrono>
#include <cstdlib>
#include <stdexcept>

#include <sys/times.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "profile_timer.h"
#include "profile_util.h"

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
  void enable()
  { _enabled = true; }

  bool enabled()
  { return _enabled; }

  void start()
  { _begin = std::chrono::high_resolution_clock::now(); }

  void stop(double *t) const
  {
    auto delta = std::chrono::high_resolution_clock::now() - _begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(delta);

    *t = static_cast<double>(ms.count()) / 10e3;
  }

private:
  bool _enabled = false;
  std::chrono::high_resolution_clock::time_point _begin;

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

pid_t timer_start() {
  if (realtime_timer.enabled())
    realtime_timer.start();
  else
    child_timer.start();

  return fork();
}

double timer_stop(pid_t child) {
  if (!stop_child(child)) {
    throw std::runtime_error("the forked child process terminated prematurely");
    return 0.0;
  }

  double t;
  if (realtime_timer.enabled())
    realtime_timer.stop(&t);
  else
    child_timer.stop(&t);

  return t;
}
