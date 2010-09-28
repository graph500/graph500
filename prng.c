/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stddef.h>
#include <stdlib.h>
#include <errno.h>

#if defined(DONT_USE_INTERNAL)
#if defined(__MTA__)
#include <mta_rng.h>

void
init_random (void)
{
  short seed[54];

  if (getenv ("SEED")) {
    errno = 0;
    seed[0] = strtol (getenv ("SEED"), NULL, 10);
    if (!errno) {
      int k;
      for (k = 1; k < 54; ++k) seed[k] = 0;
      prand_seed (seed);
    }
  }
}

void
random_vector (double *R, int64_t n)
{
  prand (n, R);
}
#else
void
init_random (void)
{
  long seed = -1;
  if (getenv ("SEED")) {
    errno = 0;
    seed = strtol (getenv ("SEED"), NULL, 10);
    if (errno) seed = -1;
  }

  if (seed < 0) seed = 0xDECAFBAD;
  srand48 (seed);
}

void
random_vector (double *R, int64_t n)
{
  int64_t k;
  OMP("omp single")
  for (k = 0; k < n; ++k)
    R[k] = drand48 ();
}
#endif
#else /* Use the generator's prng */
#include "generator/splittable_mrg.h"

uint_fast32_t prng_seed[5];
static mrg_state prng_state_store;
void *prng_state = &prng_state_store;

/* Spread the two 64-bit numbers into five nonzero values in the correct
 * range. */
static void make_mrg_seed(uint64_t userseed, uint_fast32_t* seed) {
  seed[0] = (userseed & 0x3FFFFFFF) + 1;
  seed[1] = ((userseed >> 30) & 0x3FFFFFFF) + 1;
  seed[2] = (userseed & 0x3FFFFFFF) + 1;
  seed[3] = ((userseed >> 30) & 0x3FFFFFFF) + 1;
  seed[4] = ((userseed >> 60) << 4) + (userseed >> 60) + 1;
}

void
init_random (void)
{
  long seed = -1;
  if (getenv ("SEED")) {
    errno = 0;
    seed = strtol (getenv ("SEED"), NULL, 10);
    if (errno) seed = -1;
  }

  if (seed < 0) seed = 0xDECAFBAD;
  make_mrg_seed (seed, prng_seed);
  mrg_seed(&prng_state_store, prng_seed[0], prng_seed[1], prng_seed[2],
	   prng_seed[3], prng_seed[4]);
}

void
random_vector (double *R, int64_t n)
{
  int64_t k;
  if (omp_get_num_threads () > 1) {
    OMP("omp for")
      for (k = 0; k < n; ++k) {
	mrg_state new_st = prng_state_store;
	mrg_skip(&new_st, 1, k, 0);
	R[k] = mrg_get_double_orig (&new_st);
      }
    mrg_skip(&prng_state_store, 1, n, 0);
  } else {
    for (k = 0; k < n; ++k)
      R[k] = mrg_get_double_orig (&prng_state_store);
  }
}
#endif
