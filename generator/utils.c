/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif
#include "splittable_mrg.h"
#include "graph_generator.h"
#include "permutation_gen.h"
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef __MTA__
#include <sys/mta_task.h>
#endif
#ifdef GRAPH_GENERATOR_MPI
#include <mpi.h>
#endif
#ifdef GRAPH_GENERATOR_OMP
#include <omp.h>
#endif
#include "utils.h"

void* xmalloc(size_t n) {
  void* p = malloc(n);
  if (!p) {
    fprintf(stderr, "Out of memory trying to allocate %zu byte(s)\n", n);
    abort();
  }
  return p;
}

void* xcalloc(size_t n, size_t k) {
  void* p = calloc(n, k);
  if (!p) {
    fprintf(stderr, "Out of memory trying to allocate %zu byte(s)\n", n);
    abort();
  }
  return p;
}

/* Get a number in [0, n) in an unbiased way. */
#ifdef __MTA__
#pragma mta inline
#endif
uint_fast64_t random_up_to(mrg_state* st, uint_fast64_t n) {
  /* PRNG returns values in [0, 0x7FFFFFFF) */
  /* Two iters returns values in [0, 0x3FFFFFFF00000001) */
  assert (n > 0 && n <= UINT64_C(0x3FFFFFFF00000001));
  if (n == 1) {
    return 0;
  } else if (n <= UINT64_C(0x7FFFFFFF)) {
    uint_fast64_t acc_value_limit = (UINT64_C(0x7FFFFFFF) / n) * n; /* Round down to multiple of n */
    while (1) {
      uint_fast64_t acc = mrg_get_uint_orig(st);
      if (acc >= acc_value_limit) continue;
      return acc % n;
    }
  } else if (n <= UINT64_C(0x3FFFFFFF00000001)) {
    uint_fast64_t acc_value_limit = (UINT64_C(0x3FFFFFFF00000001) / n) * n; /* Round down to multiple of n */
    while (1) {
      uint_fast64_t acc = mrg_get_uint_orig(st) * UINT64_C(0x7FFFFFFF);
      acc += mrg_get_uint_orig(st); /* Do this separately to get fixed ordering. */
      if (acc >= acc_value_limit) continue;
      return acc % n;
    }
  } else {
    /* Should have been caught before */
    return 0;
  }
}

/* Spread the two 64-bit numbers into five nonzero values in the correct
 * range. */
void make_mrg_seed(uint64_t userseed1, uint64_t userseed2, uint_fast32_t* seed) {
  seed[0] = (userseed1 & 0x3FFFFFFF) + 1;
  seed[1] = ((userseed1 >> 30) & 0x3FFFFFFF) + 1;
  seed[2] = (userseed2 & 0x3FFFFFFF) + 1;
  seed[3] = ((userseed2 >> 30) & 0x3FFFFFFF) + 1;
  seed[4] = ((userseed2 >> 60) << 4) + (userseed1 >> 60) + 1;
}


