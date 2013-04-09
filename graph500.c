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
#include "graph500-impl.h"
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

static packed_edge * restrict IJ = NULL;

static void run_bfs (void);

int
main (int argc, char **argv)
{
  get_options (argc, argv);

  init_prng ();

  if (VERBOSE)
    fprintf (stderr, "Running with %" PRId64 " vertices and %" PRId64 " edges.\n",
	     NV, NE);

  generation_time = construction_time = 0;

  /*
    If running the benchmark under an architecture simulator, replace
    the following if () {} else {} with a statement pointing IJ
    to wherever the edge list is mapped into the simulator's memory.
  */
#if defined(STORED_EDGELIST)
  size_t sz = NE * sizeof (*IJ);
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
#else
  if (dumpname) {
    fprintf (stderr, "Reading the edge list but not storing it is unsupported."
	     "  Sorry.  But not really.\n");
    return EXIT_FAILURE;
  }
#endif

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
  TIME(construction_time, err = create_graph_from_edgelist (IJ, NE, NV));
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
    sample_roots (bfs_root, NROOT, NE);
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
    int64_t *bfs_tree;

    /* Re-allocate. Some systems may randomize the addres... */
    bfs_tree = xmalloc_large (NV * sizeof (*bfs_tree));
    assert (bfs_root[m] < NV);

    if (VERBOSE) fprintf (stderr, "Running bfs %d from %" PRId64 "...", m,
			  bfs_root[m]);
    TIME(bfs_time[m], err = make_bfs_tree (bfs_tree, bfs_root[m]));
    if (VERBOSE) fprintf (stderr, " done\n");

    if (err) {
      perror ("make_bfs_tree failed");
      abort ();
    }

    if (VERBOSE) fprintf (stderr, "Verifying bfs %d...", m);
    TIME(bfs_verify_time[m],
	 bfs_depth[m] = verify_bfs_tree (bfs_tree, NULL, bfs_root[m], bfs_time[m],
					 IJ, NE));
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
