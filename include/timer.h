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

namespace timer
{

class Timer
{
  using count_type = unsigned long long;

public:
  enum Precision { SECONDS, MILLISECONDS, MICROSECONDS };

  enum { DECIMALS = 3 };

  static bool enabled;
  static std::ostream *out;

  Timer(char const *name, Precision precision = SECONDS)
  : _name(name),
    _precision(precision),
    _start(time())
  {}

  static bool exists(char const *name)
  { return _timers.find(name) != _timers.end(); }

  static void create(char const *name, Precision precision = SECONDS)
  {
    if (!exists(name))
      _timers.insert({name, Timer(name, precision)});
  }

  static void destroy(char const *name)
  { _timers.erase(find(name)); }

  static Timer *get(char const *name)
  { return &find(name)->second; }

  static std::vector<Timer *> get_all()
  {
    std::vector<Timer *> res;

    for (auto &it : _timers)
      res.push_back(&it.second);

    return res;
  }

  void start()
  { _start = time(); }

  void stop()
  {
    auto delta = time() - _start;

    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(delta);

    _meas.push_back(ns.count());
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

  std::vector<count_type>::size_type invoked() const
  { return _meas.size(); }

  double total() const
  { return scale(std::accumulate(_meas.begin(), _meas.end(), 0ULL)); }

  void mean_stddev(double *mean, double *stddev) const
  {
    count_type mean_, stddev_;

    util::mean_stddev(_meas, &mean_, &stddev_);

    *mean = scale(mean_);
    *stddev = scale(stddev_);
  }

private:
  static std::map<std::string, Timer>::iterator find(char const *name)
  {
    auto it = _timers.find(name);

    if (it == _timers.end())
      throw std::logic_error("timer does not exist");

    return it;
  }

  static std::chrono::high_resolution_clock::time_point time()
  { return std::chrono::high_resolution_clock::now(); }

  double scale(count_type ns) const
  {
    double nsd = static_cast<double>(ns);

    switch (_precision) {
      case SECONDS:
        return nsd / 1e9;
      case MILLISECONDS:
        return nsd / 1e6;
      case MICROSECONDS:
        return nsd / 1e3;
    }

    throw std::logic_error("unreachable");
  }

  char const *_name;
  Precision _precision;

  std::chrono::high_resolution_clock::time_point _start;
  std::vector<count_type> _meas;

  static std::map<std::string, Timer> _timers;
};

inline std::ostream &operator<<(std::ostream &s, Timer const &timer)
{
  s << "TIMER (" << timer.name() << "): ";

  if (timer.invoked() == 0) {
    s << " never invoked";
  } else {
    s << std::setprecision(Timer::DECIMALS);

    if (timer.invoked() == 1) {
      s << timer.total() << timer.unit();
    } else {
      double mean, stddev;
      timer.mean_stddev(&mean, &stddev);

      s << "total: " << timer.total() << timer.unit()
          << " (" << timer.invoked() << " invokations)"
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

#define TIMER_CREATE_WITH_PRECISION(name, precision) \
  TIMER_OP(timer::Timer::create(name, precision))
#define TIMER_START(name) \
  TIMER_OP(timer::Timer::create(name); timer::Timer::get(name)->start())
#define TIMER_STOP(name) \
  TIMER_OP(timer::Timer::get(name)->stop())
#define TIMER_EXISTS(name) \
  timer::Timer::exists(name)
#define TIMER_DUMP(name) \
  TIMER_OP(*timer::Timer::out << *timer::Timer::get(name) << std::endl; \
           timer::Timer::destroy(name))
#define TIMER_DUMP_ALL() \
  TIMER_OP(for (auto *t : timer::Timer::get_all()) { \
             *timer::Timer::out << *t << std::endl; \
             timer::Timer::destroy(t->name()); \
           })

#define TIMER_SECONDS timer::Timer::SECONDS
#define TIMER_MILLISECONDS timer::Timer::MILLISECONDS
#define TIMER_MICROSECONDS timer::Timer::MICROSECONDS

#endif

#endif // _GUARD_TIMER_H
