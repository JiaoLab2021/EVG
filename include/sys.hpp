#pragma once

#include <sys/resource.h>
#include <sys/time.h>

// CPU time of this process in seconds.
inline double cputime() {
    rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

// Peak resident set size of this process in bytes.
inline long peakrss() {
    rusage r;
    getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
    return r.ru_maxrss * 1024;
#else
    return r.ru_maxrss;
#endif
}

// Wall-clock time in seconds.
inline double realtime() {
    timeval tp;
    gettimeofday(&tp, nullptr);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}
