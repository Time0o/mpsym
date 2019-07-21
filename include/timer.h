#ifndef _GUARD_TIMER_H
#define _GUARD_TIMER_H

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "dbg.h"
#include "timer.h"
#include "util.h"

#ifdef NTIMER

#define Timer_create(name, precision)
#define Timer_start(name)
#define Timer_stop(name)
#define Timer_dump(name)

#else

#ifdef TIMER_CPU
#warning "Available CPU timer has low resolution"
#endif

class Timer
{
public:
  enum Precision { SECONDS, MILLISECONDS, MICROSECONDS };

  enum { DECIMALS = 3 };

  static bool enabled;
  static std::ostream &out;

  Timer(char const *name, Precision precision)
  : _start(time()),
    _name(name),
    _precision(precision)
  {}

  static void create(char const *name, Precision precision)
  {
    auto it = _timers.find(name);
    if (it != _timers.end())
      throw std::logic_error("timer already exists");

    _timers.insert({name, Timer(name, precision)});
  }

  static Timer &get(char const *name)
  {
    auto it = _timers.find(name);
    if (it == _timers.end())
      throw std::logic_error("no such timer");

    return it->second;
  }

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

  static std::unordered_map<std::string, Timer> _timers;
};

inline std::ostream &operator<< (std::ostream &s, Timer const &timer)
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

#define Timer_create(name, precision) Timer::create(name, precision)
#define Timer_start(name) Timer::get(name).start()
#define Timer_stop(name) Timer::get(name).stop()
#define Timer_dump(name) \
  do { if (Timer::enabled) Timer::out << Timer::get(name) << std::endl; } while (0)

#endif

#endif // _GUARD_TIMER_H
