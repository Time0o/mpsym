#include <atomic>
#include <chrono>
#include <cstdlib>
#include <thread>

#include "gmock/gmock.h"

#include "timeout.hpp"

#include "test_main.cpp"

using namespace mpsym::internal::timeout;

using ms = std::chrono::milliseconds;

namespace
{

template<typename REP, typename PERIOD>
void sleep(std::chrono::duration<REP, PERIOD> const &time)
{ std::this_thread::sleep_for(time); }

template<typename REP, typename PERIOD>
bool within(std::chrono::duration<REP, PERIOD> const &duration,
            std::chrono::duration<REP, PERIOD> const &target,
            std::chrono::duration<REP, PERIOD> const &margin)
{ return target - margin <= duration && target + margin >= duration; }

TEST(TimeoutTest, CanTimeoutFunction)
{
  enum : int { ARG = 42 };

  auto id_no_timeout = [](int arg)
  {
    sleep(ms(10));
    return arg;
  };

  auto id_no_timeout_timeout(ms(100));

  EXPECT_EQ(ARG, run_with_timeout("id_no_timeout", id_no_timeout_timeout, id_no_timeout, ARG))
    << "Function returns before timeout.";

  EXPECT_TRUE(within(id_no_timeout_timeout, ms(10), ms(5)))
    << "Approximate function execution time correctly determined.";

  auto id_timeout = [](int arg)
  {
    sleep(ms(1000));
    return arg;
  };

  auto id_timeout_timeout(ms(100));

  EXPECT_THAT(
    [&](){ run_with_timeout("id_timeout", id_timeout_timeout, id_timeout, ARG); },
    testing::ThrowsMessage<TimeoutError>("id_timeout timeout"))
      << "Function timeout yields exception.";

  auto endless_loop = [](std::atomic<bool> &abort, std::atomic<bool> &done)
  {
    for (;;) {
      if (abort.load())
        break;

      sleep(ms(10));
    }

    done.store(true);

    throw AbortedError("endless_loop_abort");
  };

  auto endless_loop_timeout(ms(100));

  std::atomic<bool> endless_loop_abort(false);
  std::atomic<bool> endless_loop_done(false);

  EXPECT_THAT(
    [&]() {
      run_with_timeout("endless_loop",
                       endless_loop_timeout,
                       endless_loop,
                       endless_loop_abort,
                       endless_loop_done);
    },
    testing::ThrowsMessage<TimeoutError>("endless_loop timeout"))
      << "Function timeout yields exception.";

  endless_loop_abort.store(true);

  sleep(ms(100));

  EXPECT_TRUE(endless_loop_done.load())
    << "Timed out thread terminates execution after abort flag is set.";
}

} // anonymous namespace
