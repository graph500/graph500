#define __STDC_CONSTANT_MACROS
#include <stdint.h>
#include "Random123/threefry.h"
#include "Random123/ReinterpretCtr.hpp"
#include "Random123/u01fixedpt.h"
#include <oct.h>

#define UINT32_MAX std::numeric_limits<uint32_t>::max()

/* PRNG section */
using namespace r123;
typedef ReinterpretCtr<r123array4x32, Threefry4x32> G;
static const G::key_type key = {{0xdeadbeef, 0xdecea5ed, 0x0badcafe, 0x5ca1ab1e}};

static void
iprng (int64_t v1, int64_t v2, uint32_t *ov)
{
  G::ctr_type ctr;
  G::ctr_type ctr_out;

  ctr.v[0] = int32_t (v1 >> 32);
  ctr.v[1] = int32_t (v1 & 0xFFFFFFFF);
  ctr.v[2] = int32_t (v2 >> 32);
  ctr.v[3] = int32_t (v2 & 0xFFFFFFFF);

  ctr_out = threefry4x32 (ctr, key);

  ov[0] = ctr_out.v[0];
  ov[1] = ctr_out.v[1];
  ov[2] = ctr_out.v[2];
  ov[3] = ctr_out.v[3];
}

DEFUN_DLD(PRNGCHECK, args, ,
          "PRNGCHECK")
{
  const int64_t v1 = args(0).int64_scalar_value();
  const int64_t v2 = args(1).int64_scalar_value();

  uint32_t ov[4];

  iprng (v1, v2, ov);

  return octave_value ((int32_t)ov[0]);
}

DEFUN_DLD(PRNG, args, ,
          "PRNG")
{
  const int64_t v1 = args(0).int64_scalar_value();
  const int64_t v2 = args(1).int64_scalar_value();

  FloatNDArray out (dim_vector (4, 1), -1);
  uint32_t ov[4];

  iprng (v1, v2, ov);

  out(0) = u01fixedpt_open_open_32_24 (ov[0]);
  out(1) = u01fixedpt_open_open_32_24 (ov[1]);
  out(2) = u01fixedpt_open_open_32_24 (ov[2]);
  out(3) = u01fixedpt_open_open_32_24 (ov[3]);

  return octave_value (out);
}

DEFUN_DLD(dpPRNG, args, ,
          "dpPRNG")
{
  const int64_t v1 = args(0).int64_scalar_value();
  const int64_t v2 = args(1).int64_scalar_value();

  NDArray out (dim_vector (2, 1), -1);
  uint32_t ov[4];

  iprng (v1, v2, ov);

  out(0) = u01fixedpt_open_open_64_53 ( uint64_t (ov[1]) << 32 | uint64_t (ov[0]));
  out(1) = u01fixedpt_open_open_64_53 ( uint64_t (ov[3]) << 32 | uint64_t (ov[2]));

  return octave_value (out);
}


/* Scrambling section */
#define UINT32_MAX std::numeric_limits<uint32_t>::max()

/* bitreverse and scramble taken from generator/graph_generator.c */

/* Reverse bits in a number; this should be optimized for performance
 * (including using bit- or byte-reverse intrinsics if your platform has them).
 * */
static inline uint64_t bitreverse(uint64_t x) {
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3)
#define USE_GCC_BYTESWAP /* __builtin_bswap* are in 4.3 but not 4.2 */
#endif

#ifdef FAST_64BIT_ARITHMETIC

  /* 64-bit code */
#ifdef USE_GCC_BYTESWAP
  x = __builtin_bswap64(x);
#else
  x = (x >> 32) | (x << 32);
  x = ((x >> 16) & UINT64_C(0x0000FFFF0000FFFF)) | ((x & UINT64_C(0x0000FFFF0000FFFF)) << 16);
  x = ((x >>  8) & UINT64_C(0x00FF00FF00FF00FF)) | ((x & UINT64_C(0x00FF00FF00FF00FF)) <<  8);
#endif
  x = ((x >>  4) & UINT64_C(0x0F0F0F0F0F0F0F0F)) | ((x & UINT64_C(0x0F0F0F0F0F0F0F0F)) <<  4);
  x = ((x >>  2) & UINT64_C(0x3333333333333333)) | ((x & UINT64_C(0x3333333333333333)) <<  2);
  x = ((x >>  1) & UINT64_C(0x5555555555555555)) | ((x & UINT64_C(0x5555555555555555)) <<  1);
  return x;

#else

  /* 32-bit code */
  uint32_t h = (uint32_t)(x >> 32);
  uint32_t l = (uint32_t)(x & UINT32_MAX);
#ifdef USE_GCC_BYTESWAP
  h = __builtin_bswap32(h);
  l = __builtin_bswap32(l);
#else
  h = (h >> 16) | (h << 16);
  l = (l >> 16) | (l << 16);
  h = ((h >> 8) & UINT32_C(0x00FF00FF)) | ((h & UINT32_C(0x00FF00FF)) << 8);
  l = ((l >> 8) & UINT32_C(0x00FF00FF)) | ((l & UINT32_C(0x00FF00FF)) << 8);
#endif
  h = ((h >> 4) & UINT32_C(0x0F0F0F0F)) | ((h & UINT32_C(0x0F0F0F0F)) << 4);
  l = ((l >> 4) & UINT32_C(0x0F0F0F0F)) | ((l & UINT32_C(0x0F0F0F0F)) << 4);
  h = ((h >> 2) & UINT32_C(0x33333333)) | ((h & UINT32_C(0x33333333)) << 2);
  l = ((l >> 2) & UINT32_C(0x33333333)) | ((l & UINT32_C(0x33333333)) << 2);
  h = ((h >> 1) & UINT32_C(0x55555555)) | ((h & UINT32_C(0x55555555)) << 1);
  l = ((l >> 1) & UINT32_C(0x55555555)) | ((l & UINT32_C(0x55555555)) << 1);
  return ((uint64_t)l << 32) | h; /* Swap halves */

#endif
}

/* Apply a permutation to scramble vertex numbers; a randomly generated
 * permutation is not used because applying it at scale is too expensive. */
static int64_t scramble(int64_t v0, int lgN, uint64_t val0, uint64_t val1) {
  uint64_t v = (uint64_t)v0;
  v += val0 + val1;
  v *= (val0 | UINT64_C(0x4519840211493211));
  v = (bitreverse(v) >> (64 - lgN));
  assert ((v >> lgN) == 0);
  v *= (val1 | UINT64_C(0x3050852102C843A5));
  v = (bitreverse(v) >> (64 - lgN));
  assert ((v >> lgN) == 0);
  return (int64_t)v;
}

static int scramble_already_init = 0;
static uint64_t seed0, seed1;
static void
scramble_init (void)
{
  uint32_t v[4];

  iprng (uint64_t (-1), uint64_t (-1), v);
  seed0 = uint64_t (v[0]) << 32 | uint64_t (v[1]);
  seed1 = uint64_t (v[2]) << 32 | uint64_t (v[3]);
  scramble_already_init = 1;
}

DEFUN_DLD(scramble, args, ,
          "Scramble vertex labels to reduce locality.")
{
  const int64NDArray v = args(0).int64_array_value();
  const dim_vector dims = v.dims ();
  const size_t N = v.length ();
  const int scale = args(1).int_value();

  int64NDArray out (dims, -1);

  if (!scramble_already_init) scramble_init ();

#pragma omp for
  for (size_t k = 0; k < N; ++k)
    out.elem(k) = scramble (v.elem(k), scale, seed0, seed1);

  return octave_value (out);
}
