#include <iostream>
#include <string>
#include <unordered_map>

#include "timer.h"

#ifndef NTIMER

bool Timer::enabled = false;
std::ostream &Timer::out = std::cout;
std::unordered_map<std::string, Timer> Timer::_timers;

#endif
