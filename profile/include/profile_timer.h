#ifndef _GUARD_PROFILE_TIMER_H
#define _GUARD_PROFILE_TIMER_H

#include <sys/types.h>

void timer_realtime_enable();
bool timer_realtime_enabled();

pid_t timer_start();
double timer_stop(pid_t child);

#endif // _GUARD_PROFILE_TIMER_H
