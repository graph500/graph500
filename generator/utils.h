/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef UTILS_H
#define UTILS_H

#include <stddef.h>
#include <stdint.h>
#include "splittable_mrg.h"
#ifdef __MTA__
#include <sys/mta_task.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

void* xmalloc(size_t n);
void* xcalloc(size_t n, size_t k);
uint_fast64_t random_up_to(mrg_state* st, uint_fast64_t n);
void make_mrg_seed(uint64_t userseed1, uint64_t userseed2, uint_fast32_t* seed);

/* Compare-and-swap; return 1 if successful or 0 otherwise. */
#ifdef __MTA__
#pragma mta inline
static inline int int64_t_cas(volatile int64_t* p, int64_t oldval, int64_t newval) {
  int64_t val = readfe(p);
  if (val == oldval) {
    writeef(p, newval);
    return 1;
  } else {
    writeef(p, val);
    return 0;
  }
}
#elif defined(GRAPH_GENERATOR_MPI) || defined(GRAPH_GENERATOR_SEQ)
/* Sequential */
static inline int int64_t_cas(int64_t* p, int64_t oldval, int64_t newval) {
  if (*p == oldval) {
    *p = newval;
    return 1;
  } else {
    return 0;
  }
}
#elif defined(GRAPH_GENERATOR_OMP)
/* GCC intrinsic */
static inline int int64_t_cas(volatile int64_t* p, int64_t oldval, int64_t newval) {
  return __sync_bool_compare_and_swap(p, oldval, newval);
}
#else
#error "Need to define int64_t_cas() for your system"
#endif

#ifdef __cplusplus
}
#endif

#endif /* UTILS_H */
