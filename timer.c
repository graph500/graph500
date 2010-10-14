/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#include <unistd.h>
#include <time.h>

#if defined(__MTA__)
#include <sys/mta_task.h>
long tic_ts;
#elif defined(__MacOSX__)
static AbsoluteTime tic_ts;
#elif defined(HAVE_MACH_ABSOLUTE_TIME)
#include <mach/mach_time.h>
static uint64_t tic_ts;
#else
#if defined(CLOCK_MONOTONIC)
#define TICTOC_CLOCK CLOCK_MONOTONIC
#define TICTOC_CLOCK_NAME "CLOCK_MONOTONIC"
#elif defined(CLOCK_REALTIME)
#define TICTOC_CLOCK CLOCK_REALTIME
#define TICTOC_CLOCK_NAME "CLOCK_REALTIME"
#else
#error "Failed to find a timing clock."
#endif

static struct timespec tic_ts;
#endif

void
tic (void)
{
#if defined(__MTA__)
  MTA("mta fence")
  tic_ts = mta_get_clock (0);
#elif defined(HAVE_MACH_ABSOLUTE_TIME)
  tic_ts = mach_absolute_time();  
#else
  clock_gettime (TICTOC_CLOCK, &tic_ts);
#endif
}

double
toc (void)
{
  double out;
#if defined(__MTA__)
  long ts;
  MTA("mta fence")
  ts = mta_get_clock (tic_ts);
  out = ((double)ts) * mta_clock_period ();
  /*fprintf (stderr, "%ld %g %g %g\n", ts, out, mta_clock_period(), mta_clock_freq());*/
#elif defined(HAVE_MACH_ABSOLUTE_TIME)
  uint64_t ts, nanosec;
  static mach_timebase_info_data_t info = {0,0};
  if (info.denom == 0) {
    mach_timebase_info(&info);
  }
  ts = mach_absolute_time();
  nanosec = (ts - tic_ts) * (info.numer / info.denom);
  out = 1.0e-9 * nanosec;
#else
  struct timespec ts;
  clock_gettime (TICTOC_CLOCK, &ts);
  out = (ts.tv_nsec - (double)tic_ts.tv_nsec) * 1.0e-9;
  out += (ts.tv_sec - (double)tic_ts.tv_sec);
#endif

  return out;
}
