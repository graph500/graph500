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

#if defined(_OPENMP) || defined(__MTA__)
static double take_double (double* p);
static void release_double (double* p, double val);
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
	     double * restrict R)
{
  int64_t * restrict A = A_in;
  int64_t k;

  assert (nper <= 2);

  random_vector (R, nelem);

  OMP("omp for") MTA("mta assert nodep")
    for (k = 0; k < nelem; ++k) {
      int k2;
      int64_t place;
      double Rk, Rplace;

      Rk = take_double (&R[k]);

      assert (Rk >= 0.0 && Rk <= 1.0);

      place = k + (int64_t)floor (Rk * (nelem - k));
      if (k != place) {
	assert (place > k);
	assert (place < nelem);

	Rplace = take_double (&R[place]);

	assert (Rplace >= 0.0 && Rplace <= 1.0);

	for (k2 = 0; k2 < nper; ++k2) {
	  int64_t t;
	  t = A[place*nper + k2];
	  A[place*nper + k2] = A[k*nper + k2];
	  A[k*nper + k2] = t;
	}

	release_double (&R[place], Rplace);
      }
      release_double (&R[k], Rk);
  }
}
#else
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

static void
permute_vertex_labels (int64_t * restrict IJ, int64_t nedge, int64_t max_nvtx,
		       double * restrict R, int64_t * restrict newlabel)
{
  int64_t k;

  OMP("omp for")
  for (k = 0; k < max_nvtx; ++k)
    newlabel[k] = k;

  randpermute (newlabel, max_nvtx, 1, R);

  OMP("omp for")
  for (k = 0; k < 2*nedge; ++k)
    IJ[k] = newlabel[IJ[k]];
}

static void
permute_edgelist (int64_t * restrict IJ, int64_t nedge, double *R)
{
  randpermute (IJ, nedge, 2, R);
}

#define NRAND(ne) (5 * SCALE * (ne))

void
rmat_edgelist (int64_t *IJ_in, int64_t nedge, int SCALE,
	       double A, double B, double C)
{
  int64_t * restrict IJ = IJ_in, * restrict iwork;
  double * restrict R;
  double D = 1.0 - (A + B + C);
 
  R = xmalloc_large_ext (NRAND(nedge) * sizeof (*R) + (1L<<SCALE) * sizeof (*iwork));
  iwork = (int64_t*)&R[NRAND(nedge)];

  OMP("omp parallel") {
    int64_t k;

    random_vector (R, NRAND(nedge));

    OMP("omp for") MTA("mta assert parallel") MTA("mta use 100 streams")
      for (k = 0; k < nedge; ++k)
	rmat_edge (&I(k), &J(k), SCALE, A, B, C, D, &R[NRAND(k)]);

    permute_vertex_labels (IJ, nedge, (1L<<SCALE), R, iwork);
    permute_edgelist (IJ, nedge, R);
  }

  xfree_large (R);
}

#if defined(_OPENMP)
#if 1&&(defined(__GNUC__)||defined(__INTEL_COMPILER))
/* XXX: These are not completely reliable. */
union punny {
  int64_t i64;
  double d;
};
int
double_cas(double* p, double oldval, double newval)
{
  union punny oldugh, newugh;
  oldugh.d = oldval;
  newugh.d = newval;
  return __sync_bool_compare_and_swap ((int64_t*)p, oldugh.i64, newugh.i64);
}
double
take_double (double *p)
{
  double oldval;

  do {
    __sync_synchronize ();
    oldval = *p;
  } while (!(oldval >= 0.0 && double_cas (p, oldval, -1.0)));
  return oldval;
}
void
release_double (double *p, double val)
{
  //assert (*p == -1.0);
  *p = val;
  __sync_synchronize ();
}
#else
/* XXX: These suffice for the above uses. */
double
take_double (double *p)
{
  double out;
  do {
    OMP("omp critical (TAKE)") {
      out = *p;
      if (out >= 0.0)
	*p = -1.0;
    }
  } while (out < 0.0);
  OMP("omp flush (p)");
  return out;
}
void
release_double (double *p, double val)
{
  //assert (*p == -1.0);
  OMP("omp critical (TAKE)") {
    *p = val;
    OMP("omp flush (p)");
  }
  return;
}
#endif
#elif defined(__MTA__)
int
double_cas(double* p, double oldval, double newval)
{
  int out = 0;
  double v;

  v = readfe (p);
  if (v == oldval) {
    v = newval;
    out = 1;
  }
  writeef (p, v);
  return out;
}
double
take_double (double *p)
{
  return readfe (p);
}
void
release_double (double *p, double val)
{
  writeef (p, val);
}
#else
/* double_cas isn't used sequentially. */
#endif
