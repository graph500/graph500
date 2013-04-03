/* -*- mode: C; fill-column: 70; -*- */
/* Copyright 2010-2011,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include <alloca.h> /* Portable enough... */
/* getopt should be in unistd.h */

#include "globals.h"
#include "packed_edge.h"

#define NSTAT 9
#define PRINT_STATS(lbl, israte)					\
  do {									\
    printf ("min_%s: %20.17e\n", lbl, stats[0]);			\
    printf ("firstquartile_%s: %20.17e\n", lbl, stats[1]);		\
    printf ("median_%s: %20.17e\n", lbl, stats[2]);			\
    printf ("thirdquartile_%s: %20.17e\n", lbl, stats[3]);		\
    printf ("max_%s: %20.17e\n", lbl, stats[4]);			\
    if (!israte) {							\
      printf ("mean_%s: %20.17e\n", lbl, stats[5]);			\
      printf ("stddev_%s: %20.17e\n", lbl, stats[6]);			\
    } else {								\
      printf ("harmonic_mean_%s: %20.17e\n", lbl, stats[7]);		\
      printf ("harmonic_stddev_%s: %20.17e\n", lbl, stats[8]);	\
    }									\
  } while (0)


static int
dcmp (const void *a, const void *b)
{
  const double da = *(const double*)a;
  const double db = *(const double*)b;
  if (da > db) return 1;
  if (db > da) return -1;
  if (da == db) return 0;
  fprintf (stderr, "No NaNs permitted in output.\n");
  abort ();
  return 0;
}

static void
statistics (double *out, double *data, int64_t n)
{
  long double s, mean;
  double t;
  int k;

  /* Quartiles */
  qsort (data, n, sizeof (*data), dcmp);
  out[0] = data[0];
  t = (n+1) / 4.0;
  k = (int) t;
  if (t == k)
    out[1] = data[k];
  else
    out[1] = 3*(data[k]/4.0) + data[k+1]/4.0;
  t = (n+1) / 2.0;
  k = (int) t;
  if (t == k)
    out[2] = data[k];
  else
    out[2] = data[k]/2.0 + data[k+1]/2.0;
  t = 3*((n+1) / 4.0);
  k = (int) t;
  if (t == k)
    out[3] = data[k];
  else
    out[3] = data[k]/4.0 + 3*(data[k+1]/4.0);
  out[4] = data[n-1];

  s = data[n-1];
  for (k = n-1; k > 0; --k)
    s += data[k-1];
  mean = s/n;
  out[5] = mean;
  s = data[n-1] - mean;
  s *= s;
  for (k = n-1; k > 0; --k) {
    long double tmp = data[k-1] - mean;
    s += tmp * tmp;
  }
  out[6] = sqrt (s/(n-1));

  s = (data[0]? 1.0L/data[0] : 0);
  for (k = 1; k < n; ++k)
    s += (data[k]? 1.0L/data[k] : 0);
  out[7] = n/s;
  mean = s/n;

  /*
    Nilan Norris, The Standard Errors of the Geometric and Harmonic
    Means and Their Application to Index Numbers, 1940.
    http://www.jstor.org/stable/2235723
  */
  s = (data[0]? 1.0L/data[0] : 0) - mean;
  s *= s;
  for (k = 1; k < n; ++k) {
    long double tmp = (data[k]? 1.0L/data[k] : 0) - mean;
    s += tmp * tmp;
  }
  s = (sqrt (s)/(n-1)) * out[7] * out[7];
  out[8] = s;
}

void
output_results (const char * implementation,
		const double generation_time,
		const double construction_time,
		const int64_t *root,
		const double *bfs_time, const int64_t *bfs_depth,
		const double *bfs_verify_time)
{
  int k;
  int64_t sz;
  double *tm;
  double *stats;
  char *s;

  tm = alloca (NROOT * sizeof (*tm));
  stats = alloca (NSTAT * sizeof (*stats));
  if (!tm || !stats) {
    perror ("Error allocating within final statistics calculation.");
    abort ();
  }

  s = getenv ("MACHINE");
  if (!s) s = "unknown";
  printf ("MACHINE: %s\n", s);
  s = getenv ("COMMENT");
  if (s) printf ("COMMENT: %s\n", s);
  printf ("IMPLEMENTATION: %s\n", implementation);

  printf ("SCALE: %d\nEDGEFACTOR: %d\nNROOT: %d\nMAXWEIGHT: %d\n",
	  SCALE, EF, NROOT, MAXWEIGHT);

  sz = NE * sizeof (packed_edge);
  printf ("TERASIZE: %.8e\n", sz/1.0e12);
  printf ("A: %e\nB: %e\n", (double)A, (double)B);
  printf ("K0TIME: %.8e\n", generation_time);
  printf ("K1TIME: %.8e\n", construction_time);

  for (k = 0; k < NROOT; ++k)
    tm[k] = NE / bfs_time[k];
  statistics (stats, tm, NROOT);
  printf ("K2TEPSMEAN: %.8e\n", stats[7]);
  printf ("K2TEPSSTDDEV: %.8e\n", stats[8]);

  printf ("K3TEPSMEAN: %.8e\n", (double)0);
  printf ("K3TEPSSTDDEV: %.8e\n", (double)0);

  printf ("\nroot,k2time,k2max,k2vtime,k3time,k3max,k3vtime\n");
  for (k = 0; k < NROOT; ++k)
    printf ("%" PRId64 ","
	    "%.8e,%" PRId64 ",%.8e,"
	    "%.8e,%" PRId64 ",%.8e\n",

	    root[k],
	    bfs_time[k], bfs_depth[k], bfs_verify_time[k],
	    -1.0, (int64_t)-1, -1.0);
}
