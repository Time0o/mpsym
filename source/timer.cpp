#include <iostream>
#include <map>
#include <string>

#include "timer.hpp"

#ifndef NTIMER

namespace mpsym
{

namespace internal
{

namespace timer
{

bool Timer::enabled = false;
std::ostream *Timer::out = &std::cout;
std::map<std::string, Timer> Timer::_timers;

} // namespace timer

} // namespace internal

} // namespace mpsym

#endif
