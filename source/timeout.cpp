#include "timeout.hpp"

#include <atomic>
#include <condition_variable>
#include <mutex>

namespace mpsym
{

namespace internal
{

namespace timeout
{

std::atomic<int> _timeout_thread_count(0);
std::condition_variable _timeout_thread_count_cv;
std::mutex _timeout_thread_count_mtx;

} // namespace timeout

} // namespace internal

} // namespace mpsym
