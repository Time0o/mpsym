#include <atomic>
#include <chrono>
#include <cstdlib>
#include <memory>
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

TEST(TimeoutTest, CanTimeoutFunction)
{
  enum : int { ARG = 42 };

  auto id_no_timeout = [](int arg)
  {
    sleep(ms(10));
    return arg;
  };

  EXPECT_EQ(ARG, run_with_timeout("id_no_timeout", ms(100), id_no_timeout, ARG))
    << "Function returns before timeout.";

  auto id_timeout = [](int arg)
  {
    sleep(ms(1000));
    return arg;
  };

  EXPECT_THAT(
    [&](){ run_with_timeout("id_timeout", ms(100), id_timeout, ARG); },
    testing::ThrowsMessage<TimeoutError>("id_timeout timeout"))
      << "Function timeout yields exception.";

  auto endless_loop = [](std::atomic<bool> &done, aborted_type aborted)
  {
    for (;;) {
      if (marked_aborted(aborted))
        break;

      sleep(ms(10));
    }

    done = true;

    throw AbortedError("endless_loop_abort");
  };

  std::atomic<bool> endless_loop_done(false);

  EXPECT_THAT(
    [&]() {
      run_abortable_with_timeout("endless_loop",
                                 ms(100),
                                 endless_loop,
                                 endless_loop_done);
    },
    testing::ThrowsMessage<TimeoutError>("endless_loop timeout"))
      << "Function timeout yields exception.";

  sleep(ms(100));

  EXPECT_TRUE(endless_loop_done)
    << "Timed out thread terminates execution after abort flag is set.";

  wait_for_timed_out_threads();
}

} // anonymous namespace
