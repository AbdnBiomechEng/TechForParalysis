/*
 * timeval.h    1.3 2003/01/14
 *
 * Defines gettimeofday, timeval, etc. for Win32
 *
 * By Wu Yongwei
 *
 * Modified by Arnaud Barré (2009/03/24)
 *  - Replacement of the depecrated variables _daylight, _timezone
 *    by _get_daylight and _get_timezone respectively
 *    These global variables have beed depectrated in Visual C++ 2005
 */

#ifndef _TIMEVAL_H
#define _TIMEVAL_H

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#include <winsock2.h>
#include <time.h>

#if defined(_MSC_VER) || defined(__BORLANDC__)
#define EPOCHFILETIME (116444736000000000i64)
#else
#define EPOCHFILETIME (116444736000000000LL)
#endif

struct timezone {
    long tz_minuteswest; /* minutes W of Greenwich */
    int tz_dsttime;     /* type of dst correction */
};

__inline int gettimeofday(struct timeval *tv, struct timezone *tz)
{
    FILETIME        ft;
    LARGE_INTEGER   li;
    __int64         t;
    static int      tzflag;

    if (tv)
    {
        GetSystemTimeAsFileTime(&ft);
        li.LowPart  = ft.dwLowDateTime;
        li.HighPart = ft.dwHighDateTime;
        t  = li.QuadPart;       /* In 100-nanosecond intervals */
        t -= EPOCHFILETIME;     /* Offset to the Epoch time */
        t /= 10;                /* In microseconds */
        tv->tv_sec  = (long)(t / 1000000);
        tv->tv_usec = (long)(t % 1000000);
    }

    if (tz)
    {
        if (!tzflag)
        {
            _tzset();
            tzflag++;
        }
        _get_timezone(&(tz->tz_minuteswest));
        tz->tz_minuteswest /= 60;
        _get_daylight(&(tz->tz_dsttime));
    }

    return 0;
}

#else  /* WIN32_LEAN_AND_MEAN */

#include <sys/time.h>

#endif /* _WIN32 */

#endif /* _TIMEVAL_H */
