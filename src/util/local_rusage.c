#ifndef LOCAL_RUSAGE_H
#define LOCAL_RUSAGE_H

#include <sys/time.h>
#include <sys/resource.h>

struct local_rusage {
  long ru_utimes; // user CPU time used, seconds
  long ru_utimems; // user CPU time used, microseconds
  long ru_stimes; // system CPU time used, seconds
  long ru_stimems; // system CPU time used, microseconds
  long ru_maxrss; // maximum resident set size
  long ru_ixrss; // integral shared memory size
  long ru_idrss; // integral unshared data size
  long ru_isrss; // integral unshared stack size
  long ru_minflt; // page reclaims (soft page faults)
  long ru_majflt; // page faults (hard page faults)
  long ru_nswap; // swaps
  long ru_inblock; // block input operations
  long ru_oublock; // block output operations
  long ru_msgsnd; // IPC messages sent
  long ru_msgrcv; // IPC messages received
  long ru_nsignals; // signals received
  long ru_nvcsw; // voluntary context switches
  long ru_nivcsw; // involuntary context switches
};

void fill_local_rusage_from_rusage(struct local_rusage *lrusage, struct rusage const usage)
{
  (*lrusage).ru_utimes = usage.ru_utime.tv_sec;
  (*lrusage).ru_utimems = usage.ru_utime.tv_usec;
  (*lrusage).ru_stimes = usage.ru_stime.tv_sec;
  (*lrusage).ru_stimems = usage.ru_stime.tv_usec;
  (*lrusage).ru_maxrss = usage.ru_maxrss;
  (*lrusage).ru_ixrss = usage.ru_ixrss;
  (*lrusage).ru_idrss = usage.ru_idrss;
  (*lrusage).ru_isrss = usage.ru_isrss;
  (*lrusage).ru_minflt = usage.ru_minflt;
  (*lrusage).ru_majflt = usage.ru_majflt;
  (*lrusage).ru_nswap = usage.ru_nswap;
  (*lrusage).ru_inblock = usage.ru_inblock;
  (*lrusage).ru_oublock = usage.ru_oublock;
  (*lrusage).ru_msgsnd = usage.ru_msgsnd;
  (*lrusage).ru_msgrcv = usage.ru_msgrcv;
  (*lrusage).ru_nsignals = usage.ru_nsignals;
  (*lrusage).ru_nvcsw = usage.ru_nvcsw;
  (*lrusage).ru_nivcsw = usage.ru_nivcsw;
}

int getrusage_local(struct local_rusage * lrusage)
{
  struct rusage usage;
  int const who = RUSAGE_SELF;
  int res = 0;

  res = getrusage(who, &usage);

  fill_local_rusage_from_rusage(lrusage, usage);

  return res;
}

#endif // LOCAL_RUSAGE_H
