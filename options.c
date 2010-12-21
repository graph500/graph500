/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
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
#include "options.h"

int VERBOSE = 0;
int use_RMAT = 0;

char *dumpname = NULL;
char *rootname = NULL;

double A = A_PARAM;
double B = B_PARAM;
double C = C_PARAM;
double D = 1.0 - (A_PARAM + B_PARAM + C_PARAM);

int NBFS = NBFS_max;

int64_t SCALE = default_SCALE;
int64_t edgefactor = default_edgefactor;

void
get_options (int argc, char **argv) {
  extern int opterr;
  extern int optopt;
  extern int optind;
  extern char *optarg;
  int c, err = 0;
  int nset = 0;
  int whichset = 0;

  if (getenv ("VERBOSE"))
    VERBOSE = 1;

  while ((c = getopt (argc, argv, "v?hRs:e:A:a:B:b:C:c:D:d:Vo:r:")) != -1)
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
	      "  R   : use R-MAT from SSCA2 (default: use Kronecker generator)\n"
	      "  s   : R-MAT scale (default %" PRId64 ")\n"
	      "  e   : R-MAT edge factor (default %" PRId64 ")\n"
	      "  A|a : R-MAT A (default %lg) >= 0\n"
	      "  B|b : R-MAT B (default %lg) >= 0\n"
	      "  C|c : R-MAT C (default %lg) >= 0\n"
	      "  D|d : R-MAT D (default %lg) >= 0\n"
	      "        Note: Setting 3 of A,B,C,D requires the arguments to sum to\n"
	      "        at most 1.  Otherwise, the parameters are added and normalized\n"
	      "        so that the sum is 1.\n"
	      "  V   : Enable extra (Verbose) output\n"
	      "  o   : Read the edge list from (or dump to) the named file\n"
	      "  r   : Read the BFS roots from (or dump to) the named file\n"
	      "\n"
	      "Outputs take the form of \"key: value\", with keys:\n"
	      "  SCALE\n"
	      "  edgefactor\n"
	      "  construction_time\n"
	      "  min_time\n"
	      "  firstquartile_time\n"
	      "  median_time\n"
	      "  thirdquartile_time\n"
	      "  max_time\n"
	      "  mean_time\n"
	      "  stddev_time\n"
	      "  min_nedge\n"
	      "  firstquartile_nedge\n"
	      "  median_nedge\n"
	      "  thirdquartile_nedge\n"
	      "  max_nedge\n"
	      "  mean_nedge\n"
	      "  stddev_nedge\n"
	      "  min_TEPS\n"
	      "  firstquartile_TEPS\n"
	      "  median_TEPS\n"
	      "  thirdquartile_TEPS\n"
	      "  max_TEPS\n"
	      "  harmonic_mean_TEPS\n"
	      "  harmonic_stddev_TEPS\n"
	      , default_SCALE, default_edgefactor,
	      A_PARAM, B_PARAM, C_PARAM,
	      (1.0 - (A_PARAM + B_PARAM + C_PARAM))
	      );
      exit (EXIT_SUCCESS);
      break;
    case 'V':
      VERBOSE = 1;
      break;
    case 'R':
      use_RMAT = 1;
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
      SCALE = strtol (optarg, NULL, 10);
      if (errno) {
	fprintf (stderr, "Error parsing scale %s\n", optarg);
	err = -1;
      }
      if (SCALE <= 0) {
	fprintf (stderr, "Scale must be non-negative.\n");
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
	fprintf (stderr, "Edge factor must be non-negative.\n");
	err = -1;
      }
      break;
    case 'A':
    case 'a':
      errno = 0;
      A = strtod (optarg, NULL);
      if (whichset & 1) {
	fprintf (stderr, "A already set\n");
	err = -1;
      }
      if (errno) {
	fprintf (stderr, "Error parsing A %s\n", optarg);
	err = -1;
      }
      if (A < 0) {
	fprintf (stderr, "A must be non-negative\n");
	err = -1;
      }
      whichset |= 1;
      ++nset;
      break;
    case 'B':
    case 'b':
      errno = 0;
      B = strtod (optarg, NULL);
      if (whichset & 2) {
	fprintf (stderr, "B already set\n");
	err = -1;
      }
      if (errno) {
	fprintf (stderr, "Error parsing B %s\n", optarg);
	err = -1;
      }
      if (B < 0) {
	fprintf (stderr, "B must be non-negative\n");
	err = -1;
      }
      whichset |= 2;
      ++nset;
      break;
    case 'C':
    case 'c':
      errno = 0;
      C = strtod (optarg, NULL);
      if (whichset & 4) {
	fprintf (stderr, "C already set\n");
	err = -1;
      }
      if (errno) {
	fprintf (stderr, "Error parsing C %s\n", optarg);
	err = -1;
      }
      if (C < 0) {
	fprintf (stderr, "C must be non-negative\n");
	err = -1;
      }
      whichset |= 4;
      ++nset;
      break;
    case 'D':
    case 'd':
      errno = 0;
      D = strtod (optarg, NULL);
      if (whichset & 8) {
	fprintf (stderr, "D already set\n");
	err = -1;
      }
      if (errno) {
	fprintf (stderr, "Error parsing D %s\n", optarg);
	err = -1;
      }
      if (D < 0) {
	fprintf (stderr, "D must be non-negative\n");
	err = -1;
      }
      whichset |= 8;
      ++nset;
      break;
    default:
      fprintf (stderr, "Unrecognized option\n");
      err = -1;
    }

  if (err)
    exit (EXIT_FAILURE);
  if (nset == 3) {
    switch (whichset) {
    case (4+2+1) :
      D = 1.0 - (A + B + C);
      break;
    case (8+2+1) :
      C = 1.0 - (A + B + D);
      break;
    case (8+4+1) :
      B = 1.0 - (A + C + D);
      break;
    case (8+4+2) :
      A = 1.0 - (B + C + D);
      break;
    default:
      fprintf (stderr, "Impossible combination of three bits...\n");
      abort ();
    }
    if (A < 0 || B < 0 || C < 0 || D < 0) {
      fprintf (stderr, "When setting three R-MAT parameters, all must be < 1.\n"
	       "  A = %lg\n  B = %lg\n  C = %lg\n  D = %lg\n",
	       A, B, C, D);
      exit (EXIT_FAILURE);
    }
  } else if (nset > 0) {
    long double sum = A;
    sum += B;
    sum += C;
    sum += D;
    A /= sum;
    B /= sum;
    C /= sum;
    D = 1.0 - (A + B + C);
  }
}
