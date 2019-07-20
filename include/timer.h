#ifndef _GUARD_TIMER_H
#define _GUARD_TIMER_H

#include <cassert>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "dbg.h"
#include "timer.h"
#include "util.h"

class Timer
{
public:
  enum Precision { SECONDS, MILLISECONDS, MICROSECONDS };

  enum { DECIMALS = 3 };

  static bool enabled;
  static std::ostream &out;

  Timer()
  { assert(false); }

  Timer(char const *name, Precision precision)
  : _name(name),
    _precision(precision),
    _start(std::clock())
  {}

  static void create(char const *name, Precision precision)
  {
    assert(_timers.find(name) == _timers.end());
    _timers[name] = Timer(name, precision);
  }

  static Timer &get(char const *name)
  {
    auto it = _timers.find(name);
    assert(it != timers.end());
    return it->second;
  }

  void start()
  { _start = std::clock(); }

  void stop()
  {
    auto delta = std::clock() - _start;

    double seconds =
      static_cast<double>(delta) / static_cast<double>(CLOCKS_PER_SEC);

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

  std::clock_t _start;
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

#ifdef NTIMER

#define Timer_create(name, precision)
#define Timer_start(name)
#define Timer_stop(name)
#define Timer_dump(name)

#else

#define Timer_create(name, precision) Timer::create(name, precision)
#define Timer_start(name) Timer::get(name).start()
#define Timer_stop(name) Timer::get(name).stop()
#define Timer_dump(name) \
  do { if (Timer::enabled) Timer::out << Timer::get(name) << '\n'; } while (0)

#endif

#endif // _GUARD_TIMER_H
