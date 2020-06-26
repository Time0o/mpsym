#ifndef GUARD_TIMER_H
#define GUARD_TIMER_H

#include <algorithm>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "util.hpp"

#ifdef NTIMER

#define TIMER_CREATE(name)
#define TIMER_CREATE_WITH_PRECISION(name, precision)
#define TIMER_START(name)
#define TIMER_STOP(name)
#define TIMER_DUMP(name)

#else

namespace mpsym
{

namespace internal
{

namespace timer
{

class Timer
{
  using count_type = unsigned long long;

public:
  enum Precision { SECONDS, MILLISECONDS, MICROSECONDS };

  enum { DECIMALS = 3, RESOLUTION_TEST_TICKS = 1000 };

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
    auto ns(time_delta(_start, time()));

    count_type ns_count = ns.count();

    _meas.push_back(ns_count);
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

  bool under_resolution() const
  {
    auto filter = [](count_type count){ return count > time_overhead(true); };

    decltype(_meas.size()) meas_valid =
      std::count_if(_meas.begin(), _meas.end(), filter);

    return meas_valid < _meas.size() / 2u;
  }

  std::vector<count_type>::size_type invoked() const
  { return _meas.size(); }

  double last(bool remove_overhead = true) const
  {
    auto meas_(meas(remove_overhead));

    if (meas_.empty())
      throw std::logic_error("never invoked");

    return scale(meas_.back());
  }

  double total(bool remove_overhead = true) const
  {
    auto meas_(meas(remove_overhead));

    return scale(std::accumulate(meas_.begin(), meas_.end(), 0ULL));
  }

  void mean_stddev(double *mean, double *stddev, bool remove_overhead = true) const
  {
    auto meas_(meas(remove_overhead));

    if (meas_.empty())
      throw std::logic_error("never invoked");

    count_type mean_, stddev_;

    util::mean_stddev(meas_, &mean_, &stddev_);

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

  static std::chrono::nanoseconds time_delta(
    std::chrono::high_resolution_clock::time_point start,
    std::chrono::high_resolution_clock::time_point end)
  { return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start); }

  static count_type time_overhead(bool add_stddev)
  {
    char const *t_overhead = "overhead test";

    if (!enabled)
      throw std::logic_error("timer overhead only available in macro mode");

    static count_type mean = 0ULL;
    static count_type stddev = 0ULL;

    if (mean == 0ULL) {
      for (int i = 0; i < RESOLUTION_TEST_TICKS; ++i) {
        create(t_overhead);
        get(t_overhead)->start();
        get(t_overhead)->stop();
      }

      auto *t = get(t_overhead);
      util::mean_stddev(t->_meas, &mean, &stddev);
      destroy(t_overhead);
    }

    return add_stddev ? mean + stddev : mean;
  }

  std::vector<count_type> meas(bool remove_overhead) const
  {
    if (!remove_overhead)
      return _meas;

    auto meas_overhead_removed(_meas);
    for (auto &t : meas_overhead_removed) {
      if (t <= time_overhead(false))
        t = 0;
      else
        t -= time_overhead(false);
    }

    return meas_overhead_removed;
  }

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
    s << "never invoked";
  } else if (timer.under_resolution()) {
    s << "under resolution";
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

} // namespace internal

} // namespace mpsym

#define TIMER_NS ::mpsym::internal::timer

#define TIMER_ENABLE() do { TIMER_NS :: Timer::enabled = true; } while (0)

#define TIMER_GET_OUT() TIMER_NS :: Timer::out
#define TIMER_SET_OUT(os) do { TIMER_NS :: Timer::out = os; } while (0)

#define TIMER_OP(op) do { if (TIMER_NS :: Timer::enabled) { op; } } while (0)

#define TIMER_CREATE_WITH_PRECISION(name, precision) \
  TIMER_OP(TIMER_NS :: Timer::create(name, precision))
#define TIMER_START(name) \
  TIMER_OP(TIMER_NS :: Timer::create(name); \
           TIMER_NS :: Timer::get(name)->start())
#define TIMER_STOP(name) \
  TIMER_OP(TIMER_NS :: Timer::get(name)->stop())
#define TIMER_EXISTS(name) \
  TIMER_NS :: Timer::exists(name)
#define TIMER_DUMP(name) \
  TIMER_OP(*TIMER_NS :: Timer::out << *TIMER_NS :: Timer::get(name) << std::endl; \
           TIMER_NS :: Timer::destroy(name))
#define TIMER_DUMP_ALL() \
  TIMER_OP(for (auto *t : TIMER_NS :: Timer::get_all()) { \
             *TIMER_NS :: Timer::out << *t << std::endl; \
             TIMER_NS :: Timer::destroy(t->name()); \
           })

#define TIMER_SECONDS TIMER_NS :: Timer::SECONDS
#define TIMER_MILLISECONDS TIMER_NS :: Timer::MILLISECONDS
#define TIMER_MICROSECONDS TIMER_NS :: Timer::MICROSECONDS

#endif

#endif // GUARD_TIMER_H
