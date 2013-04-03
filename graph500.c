/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010-2011,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include <alloca.h> /* Portable enough... */
#include <fcntl.h>
/* getopt should be in unistd.h */
#include <unistd.h>

#if !defined(__MTA__)
#include <getopt.h>
#endif

#include "graph500.h"
#include "globals.h"
#include "generator.h"
#include "verify.h"
#include "prng.h"
#include "timer.h"
#include "xalloc.h"
#include "options.h"
#include "output_results.h"

static int64_t bfs_root[NROOT_MAX];

static double generation_time;
static double construction_time;
static double bfs_time[NROOT_MAX];
static int64_t bfs_depth[NROOT_MAX];
static double bfs_verify_time[NROOT_MAX];

static packed_edge * restrict IJ;

static void run_bfs (void);

int
main (int argc, char **argv)
{
  ssize_t sz;

  get_options (argc, argv);

  init_prng ();

  if (VERBOSE)
    fprintf (stderr, "Running with %" PRId64 " vertices and %" PRId64 " edges.\n",
	     NV, NE);

  /*
    If running the benchmark under an architecture simulator, replace
    the following if () {} else {} with a statement pointing IJ
    to wherever the edge list is mapped into the simulator's memory.
  */
  sz = NE * sizeof (*IJ);
  IJ = xmalloc_large (sz);
  if (!dumpname) {
    if (VERBOSE) fprintf (stderr, "Generating edge list...");
    TIME(generation_time, make_graph (IJ));
    if (VERBOSE) fprintf (stderr, " done.\n");
  } else {
    int fd;
    if (VERBOSE) fprintf (stderr, "Reading edge list...");
    if ((fd = open (dumpname, O_RDONLY)) < 0) {
      perror ("Cannot open input graph file");
      return EXIT_FAILURE;
    }
    if (sz != read (fd, IJ, sz)) {
      perror ("Error reading input graph file");
      return EXIT_FAILURE;
    }
    close (fd);
    if (VERBOSE) fprintf (stderr, " done.\n");
  }

  run_bfs ();

  xfree_large (IJ);

  extern char IMPLEMENTATION[];
  output_results (IMPLEMENTATION,
		  generation_time, construction_time,
		  bfs_root,
		  bfs_time, bfs_depth, bfs_verify_time);

  return EXIT_SUCCESS;
}

void
run_bfs (void)
{
  int m, err;

  if (VERBOSE) fprintf (stderr, "Creating graph...");
  TIME(construction_time, err = create_graph_from_edgelist (IJ, NE));
  if (VERBOSE) fprintf (stderr, "done.\n");
  if (err) {
    fprintf (stderr, "Failure creating graph.\n");
    exit (EXIT_FAILURE);
  }

  /*
    If running the benchmark under an architecture simulator, replace
    the following if () {} else {} with a statement pointing bfs_root
    to wherever the BFS roots are mapped into the simulator's memory.
  */
  if (!rootname) {
    sample_roots (bfs_root);
  } else {
    int fd;
    ssize_t sz;
    if ((fd = open (rootname, O_RDONLY)) < 0) {
      perror ("Cannot open input BFS root file");
      exit (EXIT_FAILURE);
    }
    sz = NROOT * sizeof (*bfs_root);
    if (sz != read (fd, bfs_root, sz)) {
      perror ("Error reading input BFS root file");
      exit (EXIT_FAILURE);
    }
    close (fd);
  }

  for (m = 0; m < NROOT; ++m) {
    int64_t *bfs_tree, max_bfsvtx;

    /* Re-allocate. Some systems may randomize the addres... */
    bfs_tree = xmalloc_large (NV * sizeof (*bfs_tree));
    assert (bfs_root[m] < NV);

    if (VERBOSE) fprintf (stderr, "Running bfs %d from %" PRId64 "...", m,
			  bfs_root[m]);
    TIME(bfs_time[m], err = make_bfs_tree (bfs_tree, &max_bfsvtx, bfs_root[m]));
    if (VERBOSE) fprintf (stderr, "%" PRId64 " done\n", max_bfsvtx);

    if (err) {
      perror ("make_bfs_tree failed");
      abort ();
    }

    if (VERBOSE) fprintf (stderr, "Verifying bfs %d...", m);
    TIME(bfs_verify_time[m],
	 bfs_depth[m] = verify_bfs_tree (bfs_tree, max_bfsvtx, bfs_root[m], IJ, NE));
    if (VERBOSE) fprintf (stderr, "done\n");
    if (bfs_depth[m] < 0) {
      fprintf (stderr, "bfs %d from %" PRId64 " failed verification (%" PRId64 ")\n",
	       m, bfs_root[m], bfs_depth[m]);
      abort ();
    }

    xfree_large (bfs_tree);
  }

  destroy_graph ();
}

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

void
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
