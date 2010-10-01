/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <assert.h>

#include "xalloc.h"
#include "prng.h"
#include "generator/splittable_mrg.h"

#if defined(_OPENMP) || defined(__MTA__)
static int64_t take_i64 (volatile int64_t* p);
static void release_i64 (volatile int64_t* p, int64_t val);
#endif

/* Recursively divide a grid of N x N by four to a single point, (i, j).
   Choose between the four quadrants with probability a, b, c, and d.
   Create an edge between node i and j.
*/

MTA("mta expect parallel context")
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

#if defined(_OPENMP)||defined(__MTA__)
static void
randpermute (int64_t *A_in, int64_t nelem, int nper,
	     mrg_state * restrict st)
{
  int64_t * restrict A = A_in;
  int64_t k;

  assert (nper <= 2);

  OMP("omp for") MTA("mta assert nodep")
    for (k = 0; k < nelem; ++k) {
      int k2;
      int64_t place;
      double Rk, Rplace;
      int64_t Ak, Aplace;

      Ak = take_i64 (&A[k*nper]);

      assert (Ak >= 0);

      mrg_skip (st, 1, k, 0);
      place = k + (int64_t)floor (mrg_get_double_orig (st) * (nelem - k));
      if (k != place) {
	assert (place > k);
	assert (place < nelem);

	Aplace = take_i64 (&A[place*nper]);

	assert (Aplace >= 0);

	for (k2 = 1; k2 < nper; ++k2) {
	  int64_t t;
	  t = A[place*nper + k2];
	  A[place*nper + k2] = A[k*nper + k2];
	  A[k*nper + k2] = t;
	}

	{
	  int64_t t;
	  t = Aplace;
	  Aplace = Ak;
	  Ak = t;
	}

	release_i64 (&A[place*nper], Aplace);
      }
      release_i64 (&A[k*nper], Ak);
  }
}
#else
static void
randpermute (int64_t *A_in, int64_t nelem, int nper,
	     mrg_state * restrict st)
{
  int64_t * restrict A = A_in;
  int64_t k;

  assert (nper <= 2);

  for (k = 0; k < nelem; ++k) {
    int k2;
    int64_t place;

    place = k + (int64_t)floor (mrg_get_double_orig (st) * (nelem - k));

    if (k != place)
      for (k2 = 0; k2 < nper; ++k2) {
	int64_t t;
	t = A[place*nper + k2];
	A[place*nper + k2] = A[k*nper + k2];
	A[k*nper + k2] = t;
      }
  }
}
#endif

void
permute_vertex_labels (int64_t * restrict IJ, int64_t nedge, int64_t max_nvtx,
		       mrg_state * restrict st, int64_t * restrict newlabel)
{
  int64_t k;

  OMP("omp for")
  for (k = 0; k < max_nvtx; ++k)
    newlabel[k] = k;

  randpermute (newlabel, max_nvtx, 1, st);

  OMP("omp for")
  for (k = 0; k < 2*nedge; ++k)
    IJ[k] = newlabel[IJ[k]];
}

void
permute_edgelist (int64_t * restrict IJ, int64_t nedge, mrg_state *st)
{
  randpermute (IJ, nedge, 2, st);
}

#define NRAND(ne) (5 * SCALE * (ne))

void
rmat_edgelist (int64_t *IJ_in, int64_t nedge, int SCALE,
	       double A, double B, double C)
{
  int64_t * restrict IJ = IJ_in, * restrict iwork;
  double D = 1.0 - (A + B + C);
 
  iwork = xmalloc_large_ext ((1L<<SCALE) * sizeof (*iwork));

  OMP("omp parallel") {
    int64_t k;
#if !defined(__MTA__)
    double * restrict Rlocal = alloca (NRAND(1) * sizeof (*Rlocal));
#endif
    mrg_state new_st = *(mrg_state*)prng_state;

    OMP("omp for") MTA("mta assert parallel") MTA("mta use 100 streams")
      for (k = 0; k < nedge; ++k) {
	int k2;
#if defined(__MTA__)
	double * restrict Rlocal = alloca (NRAND(1) * sizeof (*Rlocal));
#endif
	mrg_skip (&new_st, 1, NRAND(1), 0);
	for (k2 = 0; k2 < NRAND(1); ++k2)
	  Rlocal[k2] = mrg_get_double_orig (&new_st);
	rmat_edge (&I(k), &J(k), SCALE, A, B, C, D, Rlocal);
      }

    OMP("omp single")
      mrg_skip (prng_state, 1, NRAND(nedge), 0);
    OMP("omp barrier");
    new_st = *(mrg_state*)prng_state;
    permute_vertex_labels (IJ, nedge, (1L<<SCALE), &new_st, iwork);
    OMP("omp single")
      mrg_skip (prng_state, 1, (1L<<SCALE), 0);
    OMP("omp barrier");
    new_st = *(mrg_state*)prng_state;
    permute_edgelist (IJ, nedge, &new_st);
  }
  mrg_skip (prng_state, 1, nedge, 0);

  xfree_large (iwork);
}

#if defined(_OPENMP)
#if defined(__GNUC__)||defined(__INTEL_COMPILER)
int64_t
take_i64 (volatile int64_t *p)
{
  int64_t oldval;

  do {
    oldval = *p;
  } while (!(oldval >= 0 && __sync_bool_compare_and_swap (p, oldval, -1)));
  return oldval;
}
void
release_i64 (volatile int64_t *p, int64_t val)
{
  assert (*p == -1);
  *p = val;
}
#else
/* XXX: These suffice for the above uses. */
int64_t
take_i64 (volatile int64_t *p)
{
  int64_t out;
  do {
    OMP("omp critical (TAKE)") {
      out = *p;
      if (out >= 0)
	*p = -1;
    }
  } while (out < 0);
  OMP("omp flush (p)");
  return out;
}
void
release_i64 (volatile int64_t *p, int64_t val)
{
  assert (*p == -1);
  OMP("omp critical (TAKE)") {
    *p = val;
    OMP("omp flush (p)");
  }
  return;
}
#endif
#elif defined(__MTA__)
int64_t
take_i64 (volatile int64_t *p)
{
  return readfe (p);
}
void
release_i64 (volatile int64_t *p, int64_t val)
{
  writeef (p, val);
}
#else
/* double_cas isn't used sequentially. */
#endif
