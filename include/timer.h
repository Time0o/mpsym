#ifndef _GUARD_TIMER_H
#define _GUARD_TIMER_H

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "dbg.h"
#include "util.h"

#ifdef NTIMER

#define TIMER_CREATE(name)
#define TIMER_CREATE_WITH_PRECISION(name, precision)
#define TIMER_START(name)
#define TIMER_STOP(name)
#define TIMER_DUMP(name)

#else

#ifdef TIMER_CPU
#warning "Available CPU timer has low resolution"
#endif

namespace timer
{

class Timer
{
public:
  enum Precision { SECONDS, MILLISECONDS, MICROSECONDS };

  enum { DECIMALS = 3 };

  static bool enabled;
  static std::ostream *out;

  Timer(char const *name, Precision precision)
  : _start(time()),
    _name(name),
    _precision(precision)
  {}

  static void create(char const *name, Precision precision = MILLISECONDS)
  {
    if (!exists(name))
      _timers.insert({name, Timer(name, precision)});
  }

  static void destroy(char const *name)
  { _timers.erase(find(name)); }

  static Timer &get(char const *name)
  { return find(name)->second; }

  void start()
  { _start = time(); }

  void stop()
  {
    double seconds;

    auto delta = time() - _start;

#ifdef TIMER_CPU
    seconds = static_cast<double>(delta) / static_cast<double>(CLOCKS_PER_SEC);
#else
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(delta);
    seconds = static_cast<double>(ns.count()) / 10e9;
#endif

    _meas.push_back(scale(seconds));
  }

  char const *name() const
  { return _name; }

  char const *unit() const
  {
    switch (_precision) {
      case SECONDS:
        return "s";
      case MILLISECONDS:
        return "ms";
      case MICROSECONDS:
        return "us";
    }

    throw std::logic_error("unreachable");
  }

  unsigned count() const
  { return static_cast<unsigned>(_meas.size()); }

  double total() const
  { return std::accumulate(_meas.begin(), _meas.end(), 0.0); }

  void mean_stddev(double *mean, double *stddev) const
  { util::mean_stddev(_meas, mean, stddev); }

private:
  static bool exists(char const *name)
  { return _timers.find(name) != _timers.end(); }

  static std::map<std::string, Timer>::iterator find(char const *name)
  {
    auto it = _timers.find(name);

    if (it == _timers.end())
      throw std::logic_error("timer does not exist");

    return it;
  }

#ifdef TIMER_CPU
  std::clock_t time()
  { return std::clock(); }

  std::clock_t _start;
#else
  std::chrono::high_resolution_clock::time_point time()
  { return std::chrono::high_resolution_clock::now(); }

  std::chrono::high_resolution_clock::time_point _start;
#endif

  double scale(double t) {
    switch (_precision) {
      case SECONDS:
        return t;
      case MILLISECONDS:
        return t * 10e3;
      case MICROSECONDS:
        return t * 10e6;
    }

    throw std::logic_error("unreachable");
  }

  char const *_name;
  Precision _precision;

  std::vector<double> _meas;

  static std::map<std::string, Timer> _timers;
};

inline std::ostream &operator<<(std::ostream &s, Timer const &timer)
{
  s << "TIMER (" << timer.name() << "): ";

  if (timer.count() == 0) {
    s << " never invoked";
  } else {
    s << std::setprecision(Timer::DECIMALS);

    if (timer.count() == 1) {
      s << timer.total() << timer.unit();
    } else {
      double mean, stddev;
      timer.mean_stddev(&mean, &stddev);

      s << "total: " << timer.total() << timer.unit()
          << " (" << timer.count() << " invokations)"
          << ", mean: " << mean << timer.unit()
          << ", stddev: " << stddev << timer.unit();
    }
  }

  return s;
}

} // namespace timer

#define TIMER_ENABLE() do { timer::Timer::enabled = true; } while (0)

#define TIMER_GET_OUT() timer::Timer::out
#define TIMER_SET_OUT(os) do { timer::Timer::out = os; } while (0)

#define TIMER_OP(op) do { if (timer::Timer::enabled) { op; } } while (0)

#define TIMER_CREATE(name) \
  TIMER_OP(timer::Timer::create(name))
#define TIMER_CREATE_WITH_PRECISION(name, precision) \
  TIMER_OP(timer::Timer::create(name, precision))
#define TIMER_START(name) \
  TIMER_OP(timer::Timer::get(name).start())
#define TIMER_STOP(name) \
  TIMER_OP(timer::Timer::get(name).stop())
#define TIMER_DUMP(name) \
  TIMER_OP(*timer::Timer::out << timer::Timer::get(name); \
           timer::Timer::destroy(name))

#define TIMER_SECONDS timer::Timer::SECONDS
#define TIMER_MILLISECONDS timer::Timer::MILLISECONDS
#define TIMER_MICROSECONDS timer::Timer::MICROSECONDS

#endif

#endif // _GUARD_TIMER_H
