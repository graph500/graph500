/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <errno.h>

#include "compat.h"

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

