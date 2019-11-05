#include <iostream>
#include <map>
#include <string>

#include "timer.h"

#ifndef NTIMER

bool Timer::enabled = false;
std::ostream &Timer::out = std::cout;
std::map<std::string, Timer> Timer::_timers;

#endif
