/* -*- mode: C; mode: folding; fill-column: 70; -*- */
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
#include "globals.h"
#include "generator.h"
#include "prng.h"
#include "timer.h"
#include "xalloc.h"
#include "options.h"

static struct packed_edge * restrict IJ;

static int64_t bfs_root[NROOT_MAX];
static double generation_time;
static double root_sample_time;

const char IMPLEMENTATION[] = "make-edgelist";

int
main (int argc, char **argv)
{
  int fd;
  ssize_t sz;

  if (sizeof (int64_t) < 8) {
    fprintf (stderr, "No 64-bit support.\n");
    return EXIT_FAILURE;
  }

  get_options (argc, argv);

  init_prng ();

  sz = NE * sizeof (*IJ);
  IJ = xmalloc_large (sz);

  if (VERBOSE)
    fprintf (stderr, "Running with %" PRId64 " vertices and %" PRId64 " edges.\n",
	     NV, NE);

  if (VERBOSE) fprintf (stderr, "Generating edge list...");
  TIME(generation_time, make_graph (IJ));
  if (VERBOSE) fprintf (stderr, " done in %g seconds.\n", generation_time);

  if (dumpname)
    fd = open (dumpname, O_WRONLY|O_CREAT|O_TRUNC, 0666);
  else
    fd = 1;

  if (fd < 0) {
    fprintf (stderr, "Cannot open output file : %s\n",
	     (dumpname? dumpname : "stdout"));
    return EXIT_FAILURE;
  }

  write (fd, IJ, NE * sizeof (*IJ));

  close (fd);

  if (rootname)
    fd = open (rootname, O_WRONLY|O_CREAT|O_TRUNC, 0666);
  else
    fd = -1;

  if (fd >= 0) {
    if (VERBOSE) fprintf (stderr, "Sample roots...");
    TIME(root_sample_time, sample_roots (bfs_root, NROOT, NE));
    if (VERBOSE) fprintf (stderr, " done in %g seconds.\n", root_sample_time);
    write (fd, bfs_root, NROOT * sizeof (*bfs_root));
    close (fd);
  }

  return EXIT_SUCCESS;
}
