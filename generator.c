#include <stdlib.h>
#include <limits.h>
#include <inttypes.h>

#include <assert.h>

#include "Random123/threefry.h"
#include "Random123/u01.h"

#include "globals.h"
#include "prng.h"
#include "packed_edge.h"
#include "generator.h"

#if defined(_OPENMP)
#define OMP(x) _Pragma(x)
#else
#define OMP(x)
#endif

#if defined(__MTA__)
#define MTA(x) _Pragma(x)
#else
#define MTA(x)
#endif

/* SCALE can be up to 40.
   Edge factor 16 takes 44 bits for indexing the list, leaving 20.
   Top 32 bits of an index will have at most 12 bits active.
*/

static inline int64_t
mult_big_mod (const int64_t k, const uint64_t hi, const uint64_t low)
{
  const uint64_t k_hi = ((uint64_t)k)>>32; /* At most 12 bits. */
  const uint64_t k_low = ((uint64_t)k)&0xFFFFFFFFul;
  uint64_t out = (k_low * low) % NE;
  if (k_hi) {
    uint64_t t = ((((k_hi * low)<<16)%NE)<<16)%NE;
    out = (out + t) % NE;
  }
  if (hi) {
    uint64_t t = ((((hi * k_low)<<16)%NE)<<16)%NE;
    out = (out + t) % NE;
    if (k_hi) {
      uint64_t t = ((((((k_hi * hi)<<32)%NE)<<16)%NE)<<16)%NE;
      out = (out + t) % NE;
    }
  }
  return out;
}

static inline int64_t
idx_to_loc (const int64_t k)
{
  return mult_big_mod (k, Z_hi, Z_low);
}

static inline int64_t
loc_to_idx (const int64_t kp)
{
  return mult_big_mod (kp, Zinv_hi, Zinv_low);
}

struct i64_pair {
  int64_t v1, v2;
};

void
edge_list (int64_t * restrict i, int64_t * restrict j, uint8_t * restrict w,
	   const int64_t ne_begin, const int64_t ne_len)
{
  assert (SCALE);

  OMP("omp parallel for") MTA("mta assert nodep")
    for (int64_t t = 0; t < ne_len; ++t) {
      const int64_t kp = ne_begin + t;
      const int64_t k = loc_to_idx (kp);
      make_edge (k, &i[t], &j[t], &w[t]);
    }
}

/* Replacable for system optimizations. */
struct i64_pair
toss_darts (const float * rnd)
{
  struct i64_pair v = { 0, 0 };
  for (int bitlvl = 0; bitlvl < SCALE; ++bitlvl) {
    float perturb = rnd[bitlvl];
    float dart = rnd[SCALE + bitlvl];
    float mu = NOISEFACT * (2 * perturb - 1);
    float Ap = A * (1 - 2 * mu / (1 - 2*B));
    float Bp = B * (1 + mu);
    v.v1 |= (dart >= Ap + Bp) << bitlvl;
    v.v2 |= ((dart >= Ap && dart < Ap + Bp) || dart >= Ap + 2*Bp) << bitlvl;
  }
  return v;
}

MTA("mta inline") void
make_edge (int64_t k, int64_t * restrict i, int64_t * restrict j,
	   uint8_t * restrict w)
{
  struct i64_pair v;
  uint8_t wgt;

  wgt = random_weight (k);
  assert (wgt > 0);
  assert (wgt <= MAXWEIGHT);

  if (k < NV) {
    /* Tree edge. */
    v.v1 = k/2;
    v.v2 = k+1;
  } else {
    /* RMAT edge */
    float rnd[2*SCALE_MAX]; /* Small, on stack. */
    random_edgevals (rnd, k);
    v = toss_darts (rnd);
  }

  v.v1 = scramble (v.v1);
  v.v2 = scramble (v.v2);
  assert (v.v1 >= 0);
  assert (v.v1 < NV);
  assert (v.v2 >= 0);
  assert (v.v2 < NV);

  *i = v.v1;
  *j = v.v2;
  *w = wgt;
}

MTA("mta inline") void
make_edge_endpoints (int64_t k, int64_t * restrict i, int64_t * restrict j)
{
  struct i64_pair v;

  if (k < NV) {
    /* Tree edge. */
    v.v1 = k/2;
    v.v2 = k+1;
  } else {
    /* RMAT edge */
    float rnd[2*SCALE_MAX]; /* Small, on stack. */
    random_edgevals (rnd, k);
    v = toss_darts (rnd);
  }

  v.v1 = scramble (v.v1);
  v.v2 = scramble (v.v2);
  assert (v.v1 >= 0);
  assert (v.v1 < NV);
  assert (v.v2 >= 0);
  assert (v.v2 < NV);

  *i = v.v1;
  *j = v.v2;
}

void
make_graph (packed_edge * result)
{
  packed_edge * restrict IJ = result;
  OMP("omp parallel for") MTA("mta assert nodep")
    for (int64_t kp = 0; kp < NE; ++kp) {
      const int64_t k = loc_to_idx (kp);
      int64_t i, j;
      uint8_t w;
      make_edge (k, &i, &j, &w);
      write_edge (&IJ[kp], i, j, w);
    }

#if !defined(NDEBUG)
  OMP("omp parallel for")
  for (int64_t kp = 0; kp < NE; ++kp) {
    int64_t i = get_v0_from_edge (&IJ[kp]);
    int64_t j = get_v1_from_edge (&IJ[kp]);
    assert (i >= 0);
    assert (j >= 0);
    assert (i < NV);
    assert (j < NV);
  }
#endif
}

void
packed_edge_list (packed_edge * result, const int64_t loc_begin,
		  const int64_t ne_begin, const int64_t ne_len)
{
  assert (SCALE);

  packed_edge * restrict IJ = result;
  OMP("omp parallel for") MTA("mta assert nodep")
    for (int64_t t = 0; t < ne_len; ++t) {
      const int64_t kp = ne_begin + t;
      const int64_t k = loc_to_idx (kp);
      int64_t i, j;
      uint8_t w;
      make_edge (k, &i, &j, &w);
      write_edge (&IJ[loc_begin+t], i, j, w);
    }
}
