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

using aborted_type = std::shared_ptr<std::atomic<bool>>;

inline aborted_type unaborted()
{ return std::make_shared<std::atomic<bool>>(false); }

inline void mark_aborted(aborted_type const &aborted)
{ return aborted->store(true); }

inline bool marked_aborted(aborted_type const &aborted)
{ return aborted->load(); }

template<typename FUNC, typename ...ARGS>
using ReturnType = decltype(std::declval<FUNC>()(std::declval<ARGS>()...));

template<typename FUNC, typename ...ARGS>
using AbortableReturnType =
  decltype(std::declval<FUNC>()(std::declval<ARGS>()...,
                                std::declval<aborted_type>()));

template<typename FUNC, typename ...ARGS>
using ReturnTypeWrapper = boost::optional<ReturnType<FUNC, ARGS...>>;

template<typename FUNC, typename ...ARGS, typename REP, typename PERIOD>
std::future<ReturnType<FUNC, ARGS...>> future_with_timeout(
  std::string const &what,
  std::chrono::duration<REP, PERIOD> const &timeout,
  FUNC &&f,
  ARGS &&...args)
{
  std::packaged_task<ReturnType<FUNC, ARGS...>()> task(
   [&]() { return f(std::forward<ARGS>(args)...); });

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

template<typename FUNC, typename ...ARGS, typename REP, typename PERIOD>
typename std::enable_if<std::is_void<ReturnType<FUNC, ARGS...>>::value>::type
run_with_timeout(std::string const &what,
                 std::chrono::duration<REP, PERIOD> const &timeout,
                 FUNC &&f,
                 ARGS &&...args)
{
  (void)future_with_timeout(
    what,
    timeout,
    [&](ARGS &&...args) {
       try {
         f(std::forward<ARGS>(args)...);
       } catch (AbortedError const &aborted) {}

       _dec_timeout_thread_count();
     },
     std::forward<ARGS>(args)...);
}

template<typename FUNC, typename ...ARGS, typename REP, typename PERIOD>
typename std::enable_if<!std::is_void<ReturnType<FUNC, ARGS...>>::value,
                        ReturnType<FUNC, ARGS...>>::type
run_with_timeout(std::string const &what,
                 std::chrono::duration<REP, PERIOD> const &timeout,
                 FUNC &&f,
                 ARGS &&...args)
{
  auto future(future_with_timeout(
    what,
    timeout,
    [&](ARGS &&...args) {
      ReturnTypeWrapper<FUNC, ARGS...> ret;

      try {
        ret = f(std::forward<ARGS>(args)...);
      } catch (AbortedError const &aborted) {}

       _dec_timeout_thread_count();

       return ret;
    },
    std::forward<ARGS>(args)...));

  return *future.get();
}

template<typename FUNC, typename ...ARGS, typename REP, typename PERIOD>
AbortableReturnType<FUNC, ARGS...>
run_abortable_with_timeout(std::string const &what,
                           std::chrono::duration<REP, PERIOD> const &timeout,
                           FUNC &&f,
                           ARGS &&...args)
{
  aborted_type aborted(unaborted());

  if (timeout <= std::chrono::duration<double>::zero())
    return f(std::forward<ARGS>(args)..., aborted);

  try {
    return run_with_timeout(what,
                            timeout,
                            std::forward<FUNC>(f),
                            std::forward<ARGS>(args)...,
                            aborted);

  } catch (TimeoutError const &) {
    mark_aborted(aborted);
    throw;
  }
}

} // namespace timeout

} // namespace internal

} // namespace mpsym

#endif // GUARD_TIMEOUT_H
