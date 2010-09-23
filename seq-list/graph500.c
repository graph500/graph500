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
#include <time.h>

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

#include <getopt.h>

#include "graph500.h"
#include "graph500-util.h"
#include "graph500-prng.h"
#include "graph500-rmat.h"

/** Verify a BFS tree, return volume or -1 if failed. */
static int64_t verify_bfs_tree (int64_t *bfs_tree, int64_t max_bfsvtx,
				int64_t root,
				const int64_t *IJ, int64_t nedge);

int VERBOSE = 0;

#define A_PARAM 0.57
#define B_PARAM 0.19
#define C_PARAM 0.19
/* Leaving D = 0.05.  D is never used explicitly.  */

static double A = A_PARAM;
static double B = B_PARAM;
static double C = C_PARAM;
static double D = 1.0 - (A_PARAM + B_PARAM + C_PARAM);

#define NBFS_max 64
static int NBFS = NBFS_max;

#define default_SCALE 14
#define default_edgefactor 83

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

  init_tictoc ();
  init_random ();

  nedge = nvtx_scale * edgefactor;
  /* Catch a few possible overflows. */
  assert (nedge >= nvtx_scale);
  assert (nedge >= edgefactor);

  IJ = xmalloc_large_ext (2 * nedge * sizeof (*IJ));

  rmat_edgelist (IJ, nedge, SCALE, A, B, C, D);

  if (getenv ("DUMPGRAPH")) {
    int k;
    FILE *g = fopen (getenv ("DUMPGRAPH"), "w");
    if (g) {
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

  xfree_large_ext ();

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

  while ((c = getopt (argc, argv, "v?hs:e:A:a:B:b:C:c:D:d:V")) != -1)
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
	      , SCALE, edgefactor, A, B, C, D
	      );
      exit (EXIT_SUCCESS);
      break;
    case 'V':
      VERBOSE = 1;
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

  xfree_large (has_adj, nvtx_scale * sizeof (*has_adj));

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

    xfree_large (bfs_tree, nvtx_scale * sizeof (*bfs_tree));
  }

  destroy_graph ();
}

static int
compute_levels (int64_t * level,
		int64_t nv, const int64_t * restrict bfs_tree, int64_t root)
{
  int err = 0;
  int64_t k;
  for (k = 0; k < nv; ++k) level[k] = -1;
  level[root] = 0;

  for (k = 0; k < nv; ++k) {
    if (level[k] >= 0) continue;
    if (!err && bfs_tree[k] >= 0 && k != root) {
      int64_t parent = k;
      int64_t nhop = 0;
      /* Run up the tree until we encounter an already-leveled vertex. */
      while (parent >= 0 && level[parent] < 0 && nhop < nv) {
	assert (parent != bfs_tree[parent]);
	parent = bfs_tree[parent];
	++nhop;
      }
      if (nhop >= nv) err = -1; /* Cycle. */
      if (parent < 0) err = -2; /* Ran off the end. */

      /* Now assign levels until we meet an already-leveled vertex */
      /* NOTE: This permits benign races if parallelized. */
      nhop += level[parent];
      parent = k;
      while (level[parent] < 0) {
	assert (nhop > 0);
	level[parent] = nhop--;
	parent = bfs_tree[parent];
      }
      assert (nhop == level[parent]);

      /* Internal check to catch mistakes in races... */
#if !defined(NDEBUG)
      nhop = 0;
      parent = k;
      int64_t lastlvl = level[k]+1;
      while (level[parent] > 0) {
	assert (lastlvl == 1 + level[parent]);
	lastlvl = level[parent];
	parent = bfs_tree[parent];
	++nhop;
      }
#endif
    }
  }
  return err;
}

int64_t
verify_bfs_tree (int64_t *bfs_tree_in, int64_t max_bfsvtx,
		 int64_t root,
		 const int64_t *IJ_in, int64_t nedge)
{
  int64_t * restrict bfs_tree = bfs_tree_in;
  const int64_t * restrict IJ = IJ_in;

  int err, nedge_traversed;
  int64_t * restrict seen_edge, * restrict level;

  const int64_t nv = max_bfsvtx+1;
  int64_t k;

  /*
    This code is horrifically contorted because many compilers
    complain about continue, return, etc. in parallel sections.
  */

  if (root > max_bfsvtx || bfs_tree[root] != root)
    return -999;

  err = 0;
  nedge_traversed = 0;
  seen_edge = xmalloc_large (2 * (nv) * sizeof (*seen_edge));
  level = &seen_edge[nv];

  for (k = 0; k < nv; ++k)
    seen_edge[k] = 0;

  err = compute_levels (level, nv, bfs_tree, root);

  if (err) goto done;

  for (k = 0; k < 2*nedge; k+=2) {
    const int64_t i = IJ[k];
    const int64_t j = IJ[k+1];
    int64_t lvldiff;

    if (i > max_bfsvtx && j <= max_bfsvtx) err = -10;
    if (j > max_bfsvtx && i <= max_bfsvtx) err = -11;
    if (err || i > max_bfsvtx /* both i & j are on the same side of max_bfsvtx */)
      continue;

    /* All neighbors must be in the tree. */
    if (bfs_tree[i] >= 0 && bfs_tree[j] < 0) err = -12;
    if (bfs_tree[j] >= 0 && bfs_tree[i] < 0) err = -13;
    if (err || bfs_tree[i] < 0 /* both i & j have the same sign */)
      continue;

    /* Both i and j are in the tree, count as a traversed edge.

       NOTE: This counts self-edges and repeated edges.  They're
       part of the input data.
    */
    ++nedge_traversed;
    /* Mark seen tree edges. */
    if (i != j) {
      if (bfs_tree[i] == j)
	seen_edge[i] = 1;
      if (bfs_tree[j] == i)
	seen_edge[j] = 1;
    }
    lvldiff = level[i] - level[j];
    /* Check that the levels differ by no more than one. */
    if (lvldiff > 1 || lvldiff < -1)
      err = -14;
  }

  if (err) goto done;

  /* Check that every BFS edge was seen and that there's only one root. */
  for (k = 0; k < nv; ++k)
    if (k != root) {
      if (bfs_tree[k] >= 0 && !seen_edge[k])
	err = -15;
      if (bfs_tree[k] == k)
	err = -16;
    }

 done:

  xfree_large (seen_edge, max_bfsvtx * sizeof (*seen_edge));
  if (err) return err;
  return nedge_traversed;
}
