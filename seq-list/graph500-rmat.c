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

#include "graph500-util.h"
#include "graph500-prng.h"

/* Recursively divide a grid of N x N by four to a single point, (i, j).
   Choose between the four quadrants with probability a, b, c, and d.
   Create an edge between node i and j.
*/

static void
rmat_edge (int64_t *iout, int64_t *jout, int SCALE,
	   double A, double B, double C, double D,
	   const double *rn)
{
  size_t rni = 0;
  int64_t i = 0, j = 0;
  int64_t bit = ((int64_t)1) << (SCALE-1);

  while (1) {
    const double r = rn[rni++];
    if (r > A) { /* outside quadrant 1 */
      if (r <= A + B) /* in quadrant 2 */
	j |= bit;
      else if (r <= A + B + C) /* in quadrant 3 */
	i |= bit;
      else { /* in quadrant 4 */
	j |= bit;
	i |= bit;
      }
    }
    if (1 == bit) break;

    /*
      Assuming R is in (0, 1), 0.95 + 0.1 * R is in (0.95, 1.05).
      So the new probabilities are *not* the old +/- 10% but
      instead the old +/- 5%.
    */
    A *= 0.95 + rn[rni++]/10;
    B *= 0.95 + rn[rni++]/10;
    C *= 0.95 + rn[rni++]/10;
    D *= 0.95 + rn[rni++]/10;
    /* Used 5 random numbers. */

    {
      const double norm = 1.0 / (A + B + C + D);
      A *= norm; B *= norm; C *= norm;
    }
    /* So long as +/- are monotonic, ensure a+b+c+d <= 1.0 */
    D = 1.0 - (A + B + C);
	
    bit >>= 1;
  }
  /* Iterates SCALE times. */
  *iout = i;
  *jout = j;
}

#define I(k) (IJ[2*(k)])
#define J(k) (IJ[2*(k)+1])

static void
randpermute (int64_t *A_in, int64_t nelem, int nper,
	     double * restrict R)
{
  int64_t * restrict A = A_in;
  int64_t k;

  assert (nper <= 2);

  random_vector (R, nelem);

  for (k = 0; k < nelem; ++k) {
    int k2;
    int64_t place;
    place = k + (int64_t)floor (R[k] * (nelem - k));

    for (k2 = 0; k2 < nper; ++k2) {
      int64_t t;
      t = A[place*nper + k2];
      A[place*nper + k2] = A[k*nper + k2];
      A[k*nper + k2] = t;
    }
  }
}


static void
permute_vertex_labels (int64_t * restrict IJ, int64_t nedge, int64_t max_nvtx,
		       double * restrict R)
{
  int64_t * restrict newlabel;
  int64_t k;

  newlabel = xmalloc_large (max_nvtx * sizeof (*newlabel));

  for (k = 0; k < max_nvtx; ++k)
    newlabel[k] = k;

  randpermute (newlabel, max_nvtx, 1, R);

  for (k = 0; k < 2*nedge; ++k)
    IJ[k] = newlabel[IJ[k]];

  xfree_large (newlabel, max_nvtx * sizeof (*newlabel));
}

static void
permute_edgelist (int64_t * restrict IJ, int64_t nedge, double *R)
{
  randpermute (IJ, nedge, 2, R);
}

#define NRAND(ne) (5 * SCALE * (ne))

void
rmat_edgelist (int64_t *IJ_in, int64_t nedge, int SCALE,
	       double A, double B, double C, double D)
{
  int64_t * restrict IJ = IJ_in;
  double * restrict R;
  int64_t k;

  R = xmalloc_large (NRAND(nedge) * sizeof (*R));
  random_vector (R, NRAND(nedge));

  for (k = 0; k < nedge; ++k)
    rmat_edge (&I(k), &J(k), SCALE, A, B, C, D, &R[NRAND(k)]);

  permute_vertex_labels (IJ, nedge, (1L<<SCALE), R);
  permute_edgelist (IJ, nedge, R);

  xfree_large (R, NRAND(nedge) * sizeof (*R));
}

