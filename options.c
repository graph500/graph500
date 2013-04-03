/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010-2011,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/* getopt should be in unistd.h */
#if HAVE_UNISTD_H
#include <unistd.h>
#else
#if !defined(__MTA__)
#include <getopt.h>
#endif
#endif

#include "graph500.h"
#include "globals.h"

int VERBOSE = 0;

char *dumpname = NULL;
char *rootname = NULL;

void
get_options (int argc, char **argv) {
  extern int opterr;
  extern int optopt;
  extern int optind;
  extern char *optarg;
  int c, err = 0;

  long scale = SCALE_DEFAULT, edgefactor = EF_DEFAULT,
    maxweight = MAXWEIGHT_DEFAULT, nroot = NROOT_DEFAULT;
  double a = A_DEFAULT, b = B_DEFAULT, noisefact = NOISEFACT_DEFAULT;

  if (getenv ("VERBOSE"))
    VERBOSE = 1;

  while ((c = getopt (argc, argv, "v?hs:e:w:A:a:B:b:N:n:Vo:r:")) != -1)
    switch (c) {
    case 'v':
      printf ("%s version %d\n", NAME, VERSION);
      exit (EXIT_SUCCESS);
      break;
    case 'h':
    case '?':
      printf ("Options:\n"
	      "  v   : version\n"
	      "  h|? : this message\n"
	      "  s   : scale (default %d)\n"
	      "  e   : edge factor (default %d)\n"
	      "  w   : maximum weight (default %d)\n"
	      "  A|a : A (default %g) >= 0\n"
	      "  B|b : B (default %g) >= 0\n"
	      "  N   : noise factor (default %g)\n"
	      "  n   : number of search roots (default %d)\n"
	      "  V   : Enable extra (Verbose) output\n"
	      "  o   : Read the edge list from (or dump to) the named file\n"
	      "  r   : Read the BFS roots from (or dump to) the named file\n"
	      "\n",
	      SCALE_DEFAULT, EF_DEFAULT, MAXWEIGHT_DEFAULT,
	      A_DEFAULT, B_DEFAULT, NOISEFACT_DEFAULT,
	      NROOT_DEFAULT);
      exit (EXIT_SUCCESS);
      break;
    case 'V':
      VERBOSE = 1;
      break;
    case 'o':
      dumpname = strdup (optarg);
      if (!dumpname) {
	fprintf (stderr, "Cannot copy dump file name.\n");
	err = 1;
      }
      break;
    case 'r':
      rootname = strdup (optarg);
      if (!rootname) {
	fprintf (stderr, "Cannot copy BFS root file name.\n");
	err = 1;
      }
      break;
    case 's':
      errno = 0;
      scale = strtol (optarg, NULL, 10);
      if (errno) {
	fprintf (stderr, "Error parsing scale %s\n", optarg);
	err = -1;
      }
      if (scale <= 0) {
	fprintf (stderr, "Scale must be positive.\n");
	err = -1;
      }
      if (scale > SCALE_MAX) {
	fprintf (stderr, "Scale cannot exceed %d.\n", SCALE_MAX);
	err = -1;
      }
      break;
    case 'e':
      errno = 0;
      edgefactor = strtol (optarg, NULL, 10);
      if (errno) {
	fprintf (stderr, "Error parsing edge factor %s\n", optarg);
	err = -1;
      }
      if (edgefactor <= 0) {
	fprintf (stderr, "Edge factor must be positive.\n");
	err = -1;
      }
      break;
    case 'w':
      errno = 0;
      maxweight = strtol (optarg, NULL, 10);
      if (errno) {
	fprintf (stderr, "Error parsing maximum weight %s\n", optarg);
	err = -1;
      }
      if (maxweight <= 0) {
	fprintf (stderr, "Maximum weight must be positive.\n");
	err = -1;
      }
      if (maxweight > UINT8_MAX) {
	fprintf (stderr, "Maximum weight must fit in eight bits.\n");
	err = -1;
      }
      break;
    case 'A':
    case 'a':
      errno = 0;
      a = strtod (optarg, NULL);
      if (errno) {
	fprintf (stderr, "Error parsing A %s\n", optarg);
	err = -1;
      }
      if (a < 0) {
	fprintf (stderr, "A must be non-negative\n");
	err = -1;
      }
      break;
    case 'B':
    case 'b':
      errno = 0;
      b = strtod (optarg, NULL);
      if (errno) {
	fprintf (stderr, "Error parsing B %s\n", optarg);
	err = -1;
      }
      if (b < 0) {
	fprintf (stderr, "B must be non-negative\n");
	err = -1;
      }
      break;
    case 'N':
      errno = 0;
      noisefact = strtod (optarg, NULL);
      if (errno) {
	fprintf (stderr, "Error parsing noise factor %s\n", optarg);
	err = -1;
      }
      if (noisefact < 0) {
	fprintf (stderr, "Noise factor must be non-negative\n");
	err = -1;
      }
      break;
    case 'n':
      errno = 0;
      nroot = strtol (optarg, NULL, 10);
      if (errno) {
	fprintf (stderr, "Error parsing scale %s\n", optarg);
	err = -1;
      }
      if (nroot <= 0) {
	fprintf (stderr, "Number of search roots must be positive.\n");
	err = -1;
      }
      if (nroot > NROOT_MAX) {
	fprintf (stderr, "Number of search roots cannot exceed %d.\n", NROOT_MAX);
	err = -1;
      }
      break;
    default:
      fprintf (stderr, "Unrecognized option\n");
      err = -1;
    }

  if (err)
    exit (EXIT_FAILURE);

  init_globals (scale, edgefactor, maxweight, nroot, a, b, noisefact);
}
