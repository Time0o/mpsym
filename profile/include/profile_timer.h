#ifndef _GUARD_PROFILE_TIMER_H
#define _GUARD_PROFILE_TIMER_H

#include <sys/types.h>

void timer_start();
double timer_stop(pid_t child = 0);

#endif // _GUARD_PROFILE_TIMER_H
