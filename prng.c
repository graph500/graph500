/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010, , Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stddef.h>
#include <stdlib.h>
#include <errno.h>

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
  for (k = 0; k < n; ++k)
    R[k] = drand48 ();
}
#endif
