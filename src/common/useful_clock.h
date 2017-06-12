#ifndef __USEFUL_CLOCK_H_
#define __USEFUL_CLOCK_H_

#include <sys/time.h> 
static long myclock()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec * 1000000) + tv.tv_usec;
}

#endif // __USEFUL_CLOCK_H_
