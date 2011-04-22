/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
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
#include <unistd.h>

#if !defined(__MTA__)
#include <getopt.h>
#endif


#include "graph500.h"
#include "rmat.h"
#include "kronecker.h"
#include "verify.h"
#include "prng.h"
#include "xalloc.h"
#include "options.h"
#include "generator/splittable_mrg.h"
#include "generator/make_graph.h"

static int64_t nvtx_scale;

static struct packed_edge * restrict IJ;
static int64_t nedge;

static int64_t bfs_root[NBFS_max];

int
main (int argc, char **argv)
{
  int * restrict has_adj;
  int fd;
  int64_t desired_nedge;
  if (sizeof (int64_t) < 8) {
    fprintf (stderr, "No 64-bit support.\n");
    return EXIT_FAILURE;
  }

  if (argc > 1)
    get_options (argc, argv);

  nvtx_scale = 1L<<SCALE;

  init_random ();

  desired_nedge = nvtx_scale * edgefactor;
  /* Catch a few possible overflows. */
  assert (desired_nedge >= nvtx_scale);
  assert (desired_nedge >= edgefactor);


  if (VERBOSE) fprintf (stderr, "Generating edge list...");
  if (use_RMAT) {
    nedge = desired_nedge;
    IJ = xmalloc_large_ext (nedge * sizeof (*IJ));
    rmat_edgelist (IJ, nedge, SCALE, A, B, C);
  } else {
    make_graph(SCALE, desired_nedge, userseed, userseed, &nedge, (struct packed_edge**)(&IJ));
  }
  if (VERBOSE) fprintf (stderr, " done.\n");

  if (dumpname)
    fd = open (dumpname, O_WRONLY|O_CREAT|O_TRUNC, 0666);
  else
    fd = 1;

  if (fd < 0) {
    fprintf (stderr, "Cannot open output file : %s\n",
	     (dumpname? dumpname : "stdout"));
    return EXIT_FAILURE;
  }

  write (fd, IJ, 2 * nedge * sizeof (*IJ));

  close (fd);

  if (rootname)
    fd = open (rootname, O_WRONLY|O_CREAT|O_TRUNC, 0666);
  else
    fd = -1;

  if (rootname >= 0) {
    has_adj = xmalloc_large (nvtx_scale * sizeof (*has_adj));
    OMP("omp parallel") {
      OMP("omp for")
	for (int64_t k = 0; k < nvtx_scale; ++k)
	  has_adj[k] = 0;
      MTA("mta assert nodep") OMP("omp for")
	for (int64_t k = 0; k < nedge; ++k) {
	  const int64_t i = get_v0_from_edge(&IJ[k]);
	  const int64_t j = get_v1_from_edge(&IJ[k]);
	  if (i != j)
	    has_adj[i] = has_adj[j] = 1;
	}
    }

    /* Sample from {0, ..., nvtx_scale-1} without replacement. */
    {
      int m = 0;
      int64_t t = 0;
      while (m < NBFS && t < nvtx_scale) {
	double R = mrg_get_double_orig (prng_state);
	if (!has_adj[t] || (nvtx_scale - t)*R > NBFS - m) ++t;
	else bfs_root[m++] = t++;
      }
      if (t >= nvtx_scale && m < NBFS) {
	if (m > 0) {
	  fprintf (stderr, "Cannot find %d sample roots of non-self degree > 0, using %d.\n",
		   NBFS, m);
	  NBFS = m;
	} else {
	  fprintf (stderr, "Cannot find any sample roots of non-self degree > 0.\n");
	  exit (EXIT_FAILURE);
	}
      }
    }

    xfree_large (has_adj);
    write (fd, bfs_root, NBFS * sizeof (*bfs_root));
    close (fd);
  }

  return EXIT_SUCCESS;
}
