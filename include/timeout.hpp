#ifndef GUARD_TIMEOUT_H
#define GUARD_TIMEOUT_H

#include <chrono>
#include <functional>
#include <future>
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

template<typename FUNC, typename ...ARGS>
using ReturnType = decltype(std::declval<FUNC>()(std::declval<ARGS>()...));

template<typename FUNC, typename ...ARGS>
using ReturnTypeWrapper = boost::optional<ReturnType<FUNC, ARGS...>>;

template<typename FUNC, typename ...ARGS, typename REP, typename PERIOD>
std::future<ReturnType<FUNC, ARGS...>> future_with_timeout(
  std::string const &what,
  std::chrono::duration<REP, PERIOD> &timeout,
  FUNC &&f,
  ARGS &&...args)
{
  using namespace std::chrono;

  std::packaged_task<ReturnType<FUNC, ARGS...>()> task(
   [&]() { return f(std::forward<ARGS>(args)...); });

  auto future(task.get_future());

  auto start_thread(high_resolution_clock::now());

  std::thread thread(std::move(task));

  if (future.wait_for(timeout) == std::future_status::timeout) {
     thread.detach();

     throw TimeoutError(what);

  } else {
     thread.join();

     auto end_thread(high_resolution_clock::now());

     timeout = duration_cast<duration<REP, PERIOD>>(end_thread - start_thread);

     return future;
  }
}

template<typename FUNC, typename ...ARGS, typename REP, typename PERIOD>
typename std::enable_if<std::is_void<ReturnType<FUNC, ARGS...>>::value>::type
run_with_timeout(std::string const &what,
                 std::chrono::duration<REP, PERIOD> &timeout,
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
     },
     std::forward<ARGS>(args)...);
}

template<typename FUNC, typename ...ARGS, typename REP, typename PERIOD>
typename std::enable_if<!std::is_void<ReturnType<FUNC, ARGS...>>::value,
                        ReturnType<FUNC, ARGS...>>::type
run_with_timeout(std::string const &what,
                 std::chrono::duration<REP, PERIOD> &timeout,
                 FUNC &&f,
                 ARGS &&...args)
{
  auto future(future_with_timeout(
    what,
    timeout,
    [&](ARGS &&...args) {
      try {
        return ReturnTypeWrapper<FUNC, ARGS...>(
          f(std::forward<ARGS>(args)...));
      } catch (AbortedError const &aborted) {
        return ReturnTypeWrapper<FUNC, ARGS...>();
      }
    },
    std::forward<ARGS>(args)...));

  return *future.get();
}

} // namespace timeout

} // namespace internal

} // namespace mpsym

#endif // GUARD_TIMEOUT_H
