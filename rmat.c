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

#include "compat.h"
#include "xalloc.h"
#include "prng.h"

#if defined(_OPENMP) || defined(__MTA__)
static int double_cas(double* p, double oldval, double newval);
#endif

/* Recursively divide a grid of N x N by four to a single point, (i, j).
   Choose between the four quadrants with probability a, b, c, and d.
   Create an edge between node i and j.
*/

MTA("mta expect serial context")
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

#if defined(_OPENMP)
static void
randpermute (int64_t *A_in, int64_t nelem, int nper,
	     double * restrict R)
{
  int64_t * restrict A = A_in;
  int64_t k;

  assert (nper <= 2);

  random_vector (R, nelem);

  OMP("omp for")
    for (k = 0; k < nelem; ++k) {
      int k2;
      int64_t place;
      double Rk, Rplace;

      do {
	OMP("omp flush (R)");
	Rk = R[k];
      } while (Rk < 0.0 || !double_cas (&R[k], Rk, -1.0));

      assert (R[k] == -1.0);
      assert (Rk >= 0.0 && Rk <= 1.0);

      place = k + (int64_t)floor (Rk * (nelem - k));
      if (k != place) {
	assert (place > k);
	assert (place < nelem);

	do {
	  OMP("omp flush (R)");
	  Rplace = R[place];
	} while (Rplace < 0.0 || !double_cas (&R[place], Rplace, -1.0));

	assert (R[place] == -1.0);
	assert (Rplace >= 0.0 && Rplace <= 1.0);

	for (k2 = 0; k2 < nper; ++k2) {
	  int64_t t;
	  t = A[place*nper + k2];
	  A[place*nper + k2] = A[k*nper + k2];
	  A[k*nper + k2] = t;
	}

	R[place] = Rplace;
      }
      R[k] = Rk;
      OMP("omp flush (R)");
  }
}
#elif defined(__MTA__)
static void
randpermute (int64_t *A_in, int64_t nelem, int nper,
	     double * restrict R)
{
  int64_t * restrict A = A_in;
  int64_t k;
  double Rk, Rplace;

  assert (nper <= 2);

  random_vector (R, nelem);

  MTA("mta assert nodep")
  for (k = 0; k < nelem; ++k) {
    int k2;
    int64_t place;
    double Rk;

    Rk = readfe (&R[k]);

    place = k + (int64_t)floor (R[k] * (nelem - k));

    if (k != place) {
      double Rplace;
      Rplace = readfe (&R[place]);

      for (k2 = 0; k2 < nper; ++k2) {
	int64_t t;
	t = A[place*nper + k2];
	A[place*nper + k2] = A[k*nper + k2];
	A[k*nper + k2] = t;
      }
      writeef (&R[place], Rplace);
    }

    writeef (&R[k], Rk);
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
 
  R = xmalloc_large (NRAND(nedge) * sizeof (*R) + (1L<<SCALE) * sizeof (*iwork));
  iwork = (int64_t*)&R[NRAND(nedge)];

  OMP("omp parallel") {
    int64_t k;

    random_vector (R, NRAND(nedge));

    OMP("omp for")
      for (k = 0; k < nedge; ++k)
	rmat_edge (&I(k), &J(k), SCALE, A, B, C, D, &R[NRAND(k)]);

    permute_vertex_labels (IJ, nedge, (1L<<SCALE), R, iwork);
    permute_edgelist (IJ, nedge, R);
  }

  xfree_large (R);
}

#if defined(_OPENMP)
#if defined(__GNUC__)
int
double_cas(double* p, double oldval, double newval)
{
  union ugh {
    int64_t i64;
    double d;
  } oldugh, newugh;
  oldugh.d = oldval;
  newugh.d = newval;
  return __sync_bool_compare_and_swap ((int64_t*)p, oldugh.i64, newugh.i64);
}
#else
/* XXX: These are not correct, but suffice for the above uses. */
int
double_cas(double* p, double oldval, double newval)
{
  int out = 0;
  OMP("omp critical (CAS)") {
    double v = *p;
    if (v == oldval) {
      *p = newval;
      out = 1;
    }
  }
  OMP("omp flush (p)");
  return out;
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
#else
/* double_cas isn't used sequentially. */
#endif
