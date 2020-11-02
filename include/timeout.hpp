#ifndef GUARD_TIMEOUT_H
#define GUARD_TIMEOUT_H

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>

#include <boost/optional.hpp>

namespace mpsym
{

namespace internal
{

namespace timeout
{

// keep track of and wait for termination of detached threads

extern std::atomic<int> _timeout_thread_count;
extern std::condition_variable _timeout_thread_count_cv;
extern std::mutex _timeout_thread_count_mtx;

inline void _inc_timeout_thread_count()
{ ++_timeout_thread_count; }

inline void _dec_timeout_thread_count()
{
  --_timeout_thread_count;
  _timeout_thread_count_cv.notify_one();
}

inline void wait_for_timed_out_threads()
{
  std::unique_lock<std::mutex> lock(_timeout_thread_count_mtx);

  _timeout_thread_count_cv.wait(lock,
                                [&]{ return _timeout_thread_count == 0; });
}

// exceptions

struct TimeoutError : public std::runtime_error
{
  TimeoutError(std::string const &what)
  : std::runtime_error(what + " timeout")
  {}
};

struct AbortedError : public std::runtime_error
{
  AbortedError(std::string const &what)
  : std::runtime_error(what + " aborted")
  {}
};

// timeout function wrappers

using flag = std::shared_ptr<std::atomic<bool>>;

inline flag unset()
{ return std::make_shared<std::atomic<bool>>(false); }

inline void set(flag const &f)
{ f->store(true); }

inline bool is_set(flag const &f)
{ return f->load(); }

template<typename FUNC>
using ReturnType = decltype(std::declval<FUNC>()());

template<typename FUNC>
using AbortableReturnType =
  decltype(std::declval<FUNC>()(std::declval<flag>()));

template<typename FUNC>
using ReturnTypeWrapper = boost::optional<ReturnType<FUNC>>;

template<typename FUNC, typename REP, typename PERIOD>
std::future<ReturnType<FUNC>> future_with_timeout(
  std::string const &what,
  std::chrono::duration<REP, PERIOD> const &timeout,
  FUNC &&f)
{
  std::packaged_task<ReturnType<FUNC>()> task(f);

  auto future(task.get_future());

  _inc_timeout_thread_count();

  std::thread thread(std::move(task));

  if (future.wait_for(timeout) == std::future_status::timeout) {
    thread.detach();

    throw TimeoutError(what);

  } else {
    thread.join();

    return future;
  }
}

template<typename FUNC, typename REP, typename PERIOD>
typename std::enable_if<std::is_void<ReturnType<FUNC>>::value>::type
run_with_timeout(std::string const &what,
                 std::chrono::duration<REP, PERIOD> const &timeout,
                 FUNC &&f)
{
  (void)future_with_timeout(
    what,
    timeout,
    [&]() {
       try {
         f();
       } catch (AbortedError const &aborted) {}

       _dec_timeout_thread_count();
     });
}

template<typename FUNC, typename REP, typename PERIOD>
typename std::enable_if<!std::is_void<ReturnType<FUNC>>::value,
                        ReturnType<FUNC>>::type
run_with_timeout(std::string const &what,
                 std::chrono::duration<REP, PERIOD> const &timeout,
                 FUNC &&f)
{
  auto future(future_with_timeout(
    what,
    timeout,
    [&]() {
      ReturnTypeWrapper<FUNC> ret;

      try {
        ret = f();
      } catch (AbortedError const &) {}

       _dec_timeout_thread_count();

       return ret;
    }));

  return *future.get();
}

template<typename FUNC, typename REP, typename PERIOD>
AbortableReturnType<FUNC>
run_abortable_with_timeout(std::string const &what,
                           std::chrono::duration<REP, PERIOD> const &timeout,
                           FUNC &&f)
{
  flag aborted(unset());

  if (timeout <= std::chrono::duration<double>::zero())
    return f(aborted);

  try {
    return run_with_timeout(what,
                            timeout,
                            [&]{ return f(aborted); });

  } catch (TimeoutError const &) {
    set(aborted);
    throw;
  }
}

} // namespace timeout

} // namespace internal

} // namespace mpsym

#endif // GUARD_TIMEOUT_H
