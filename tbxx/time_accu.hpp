#ifndef TBXX_TIME_ACCU_HPP
#define TBXX_TIME_ACCU_HPP

#ifdef _MSC_VER
#include <winsock2.h>

inline void timeradd(timeval *a, timeval *b, timeval *res) {
  res->tv_sec = a->tv_sec + b->tv_sec;
  res->tv_usec = a->tv_usec + b->tv_usec;
  if(res->tv_usec > 1000000) {
    ++res->tv_sec;
    res->tv_usec -= 1000000;
  }
}

inline void timersub(timeval *a, timeval *b, timeval *res) {
  res->tv_sec = a->tv_sec - b->tv_sec;
  res->tv_usec = a->tv_usec - b->tv_usec;
  if(res->tv_usec < 0) {
    --res->tv_sec;
    res->tv_usec += 1000000;
  }
}

/* C.f. the following pieces of official Microsoft documentation:

GetSystemTimeAsFileTime function:
https://msdn.microsoft.com/en-us/library/windows/desktop/ms724397(v=vs.85).aspx

FILETIME structure:
https://msdn.microsoft.com/en-us/library/windows/desktop/ms724284(v=vs.85).aspx

Conversion from Windows begining of time to Unix epoch:
http://stackoverflow.com/a/6161842/1428153

Note: time zone not supported but we don't need it for the purpose
of time_accu below.
*/
inline int gettimeofday(struct timeval *tv, void *tz) {
  union {
    // time since 1 Jan 1601 in 10-th of micro-seconds (usec)
    unsigned long long t;
    FILETIME ft;
  };
  GetSystemTimeAsFileTime (&ft);
  t /= 10LL; //in usec
  tv->tv_usec = (long)(t % 1000000);
  tv->tv_sec  = (long)(t / 1000000 - 11644473600LL);
  return 0;
}

#else

#include <sys/time.h>

#endif

namespace tbxx {

  struct time_accu
  {
    timeval accu;
    timeval mark;

    time_accu()
    {
      timerclear(&accu);
      timerclear(&mark);
    }

    void
    set_mark()
    {
      gettimeofday(&mark, 0);
    }

    time_accu&
    accumulate()
    {
      timeval curr;
      gettimeofday(&curr, 0);
      timeval diff;
      timersub(&curr, &mark, &diff);
      timerclear(&mark);
      timeval sum;
      timeradd(&accu, &diff, &sum);
      accu.tv_sec = sum.tv_sec;
      accu.tv_usec = sum.tv_usec;
      return *this;
    }

    double
    as_double() const
    {
      return static_cast<double>(accu.tv_sec)
           + static_cast<double>(accu.tv_usec) * 1.e-6;
    }
  };

} // namespace tbxx

#endif // GUARD
