/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _FILE_OFFSET_BITS 64
#define _THREAD_SAFE
#define _XOPEN_SOURCE 600
#define _XOPEN_SOURCE_EXTENDED
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#if __STDC_VERSION__ >= 199901L
#include <inttypes.h>
#else
#warning "Defining long as int64_t."
typedef long int64_t;
#define PRId64 "ld"
#define SCNd64 "ld"
#if !defined(restrict)
#define restrict
#endif
#endif

#include <alloca.h> /* Portable enough... */
#if !defined(__MTA__)
#include <getopt.h>
#endif

#include "graph500.h"
#include "rmat.h"
#include "verify.h"
#include "prng.h"
#include "timer.h"
#include "xalloc.h"

int VERBOSE = 0;
int use_RMAT = 1; /* XXX: Fix to link in Kronecker generator... */

#define A_PARAM 0.57
#define B_PARAM 0.19
#define C_PARAM 0.19
/* Hence D = 0.05. */

static double A = A_PARAM;
static double B = B_PARAM;
static double C = C_PARAM;
static double D = 1.0 - (A_PARAM + B_PARAM + C_PARAM);

#define NBFS_max 64
static int NBFS = NBFS_max;

#define default_SCALE ((int64_t)14)
#define default_edgefactor ((int64_t)83)

static int64_t SCALE = default_SCALE;
static int64_t nvtx_scale;
static int64_t edgefactor = default_edgefactor;

static int64_t bfs_root[NBFS_max];

static double construction_time;
static double bfs_time[NBFS_max];
static int64_t bfs_nedge[NBFS_max];

static int64_t * restrict IJ;
static int64_t nedge;

static void get_options (int argc, char **argv);
static void run_bfs (void);
static void output_results (const int64_t SCALE, int64_t nvtx_scale,
			    int64_t edgefactor,
			    const double A, const double B,
			    const double C, const double D,
			    const double construction_time,
			    const int NBFS,
			    const double *bfs_time, const int64_t *bfs_nedge);

int
main (int argc, char **argv)
{
  if (sizeof (int64_t) < 8) {
    fprintf (stderr, "No 64-bit support.\n");
    return EXIT_FAILURE;
  }

  if (argc > 1)
    get_options (argc, argv);

  nvtx_scale = 1L<<SCALE;

  init_random ();

  nedge = nvtx_scale * edgefactor;
  /* Catch a few possible overflows. */
  assert (nedge >= nvtx_scale);
  assert (nedge >= edgefactor);

  IJ = xmalloc_large_ext (2 * nedge * sizeof (*IJ));

  rmat_edgelist (IJ, nedge, SCALE, A, B, C);

  if (getenv ("DUMPGRAPH")) {
    int k;
    FILE *g = fopen (getenv ("DUMPGRAPH"), "w");
    if (g) {
      fprintf (g, "%" PRId64 "\n", nedge);
      for (k = 0; k < 2*nedge; k+=2) {
	const int64_t i = IJ[k];
	const int64_t j = IJ[k+1];
	fprintf (g, "%" PRId64 " %" PRId64 "\n", i, j);
      }
      fclose (g);
    } else
      fprintf (stderr, "Failed to open \"%s\" for dumping.\n", getenv ("DUMPGRAPH"));
  }

  run_bfs ();

  xfree_large (IJ);

  output_results (SCALE, nvtx_scale, edgefactor, A, B, C, D,
		  construction_time, NBFS, bfs_time, bfs_nedge);

  return EXIT_SUCCESS;
}

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

  while ((c = getopt (argc, argv, "v?hRs:e:A:a:B:b:C:c:D:d:V")) != -1)
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

void
run_bfs (void)
{
  int *has_adj;
  int k, m, err;
  int64_t t;
  double R[2*NBFS];

  TIME(construction_time, err = create_graph_from_edgelist (IJ, nedge));
  if (err) {
    fprintf (stderr, "Failure creating graph.\n");
    exit (EXIT_FAILURE);
  }

  has_adj = xmalloc_large (nvtx_scale * sizeof (*has_adj));
  for (k = 0; k < nvtx_scale; ++k)
    has_adj[k] = 0;
  for (k = 0; k < 2*nedge; k+=2) {
    const int64_t i = IJ[k];
    const int64_t j = IJ[k+1];
    if (i != j)
      has_adj[i] = has_adj[j] = 1;
  }

  /* Sample from {0, ..., nvtx_scale-1} without replacement. */
  m = 0;
  t = 0;
  while (m < NBFS) {
    random_vector (R, 2*NBFS);
    k = 0;
    while (m < NBFS && k < 2*NBFS && t < nvtx_scale) {
      if (!has_adj[t] || (nvtx_scale - t)*R[k++] > NBFS - m) ++t;
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

  for (k = 0; k < NBFS; ++k) {
    int64_t *bfs_tree, max_bfsvtx;

    /* Re-allocate. Some systems may randomize the addres... */
    bfs_tree = xmalloc_large (nvtx_scale * sizeof (*bfs_tree));
    assert (bfs_root[k] < nvtx_scale);

    TIME(bfs_time[k], err = make_bfs_tree (bfs_tree, &max_bfsvtx, bfs_root[k]));

    if (err) {
      perror ("make_bfs_tree failed");
      abort ();
    }

    bfs_nedge[k] = verify_bfs_tree (bfs_tree, max_bfsvtx, bfs_root[k], IJ, nedge);
    if (bfs_nedge[k] < 0) {
      fprintf (stderr, "bfs %d from %" PRId64 " failed verification (%" PRId64 ")\n",
	       k, bfs_root[k], bfs_nedge[k]);
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

void
output_results (const int64_t SCALE, int64_t nvtx_scale, int64_t edgefactor,
		const double A, const double B, const double C, const double D,
		const double construction_time,
		const int NBFS, const double *bfs_time, const int64_t *bfs_nedge)
{
  int k;
  int64_t sz;
  double *tm;
  double *stats;

  tm = alloca (NBFS * sizeof (*tm));
  stats = alloca (NSTAT * sizeof (*stats));
  if (!tm || !stats) {
    perror ("Error allocating within final statistics calculation.");
    abort ();
  }

  sz = (1L << SCALE) * edgefactor * sizeof (int64_t);
  printf ("SCALE: %" PRId64 "\nnvtx: %" PRId64 "\nedgefactor: %" PRId64 "\n"
	  "tebisize: %20.17e\n",
	  SCALE, nvtx_scale, edgefactor, sz/1.0e12);
  printf ("A: %20.17e\nB: %20.17e\nC: %20.17e\nD: %20.17e\n", A, B, C, D);
  printf ("construction_time: %20.17e\n", construction_time);
  printf ("nbfs: %d\n", NBFS);

  memcpy (tm, bfs_time, NBFS*sizeof(tm[0]));
  statistics (stats, tm, NBFS);
  PRINT_STATS("time", 0);

  for (k = 0; k < NBFS; ++k)
    tm[k] = bfs_nedge[k];
  statistics (stats, tm, NBFS);
  PRINT_STATS("nedge", 0);

  for (k = 0; k < NBFS; ++k)
    tm[k] = bfs_nedge[k] / bfs_time[k];
  statistics (stats, tm, NBFS);
  PRINT_STATS("TEPS", 1);
}
