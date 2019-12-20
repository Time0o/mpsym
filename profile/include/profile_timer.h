#ifndef _GUARD_PROFILE_TIMER_H
#define _GUARD_PROFILE_TIMER_H

#include <sys/types.h>

void timer_realtime_enable();
bool timer_realtime_enabled();

void timer_start();
double timer_stop(pid_t child = 0);

void debug_timer_dump(char const *timer);

#endif // _GUARD_PROFILE_TIMER_H
