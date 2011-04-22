/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stddef.h>
#include <stdlib.h>
#include <errno.h>

#include "generator/splittable_mrg.h"

uint64_t userseed;
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
  userseed = seed;
  make_mrg_seed (seed, prng_seed);
  mrg_seed(&prng_state_store, prng_seed);
}
