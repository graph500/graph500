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

#if defined(_OPENMP)
#define OMP(x_) _Pragma(x_)
#else
#define OMP(x_)
#endif

#if defined(HAVE_LIBNUMA)
#include <numa.h>
static int numa_inited = 0;
static int numa_avail = -1;

void *
xmalloc (size_t sz)
{
  void * out;
  if (!numa_inited) {
    OMP("omp critical") {
      numa_inited = 1;
      numa_avail = numa_available ();
    }
  }

  if (numa_avail)
    out = numa_alloc (sz);
  else
    out = malloc (sz);
  if (!out) {
    fprintf(stderr, "Out of memory trying to allocate %zu byte(s)\n", sz);
    abort ();
  }
  return out;
}
void *
xcalloc (size_t n, size_t sz)
{
  void * out;
  if (!numa_inited) {
    OMP("omp critical") {
      numa_inited = 1;
      numa_avail = numa_available ();
    }
  }

  if (numa_avail) {
    size_t to_alloc;
    to_alloc = n * sz;
    if (to_alloc < n || to_alloc < sz) {
      fprintf(stderr,
	      "Allocation size out of range for %zu items of %zu byte(s)\n",
	      n, sz);
      abort ();
    }
    out = numa_alloc (n * sz);
#if defined(_OPENMP)
#pragma omp parallel for
      for (size_t k = 0; k < n; ++k)
	memset (out + k * sz, 0, sz);
#else
    memset (out, 0, n * sz);
#endif
  } else
    out = calloc (n, sz);
  if (!out) {
    fprintf(stderr,
	    "Out of memory trying to allocate/clear %zu items of %zu byte(s)\n",
	    n, sz);
    abort ();
  }
  return out;
}
void
xfree (void * p, size_t sz)
{
  if (!p) return;
  if (numa_avail >= 0)
    numa_free (p, sz);
  else
    free (p);
}
#else
void *
xmalloc (size_t sz)
{
  void * out;
  out = malloc (sz);
  if (!out) {
    fprintf(stderr, "Out of memory trying to allocate %zu byte(s)\n", sz);
    abort ();
  }
  return out;
}
void *
xcalloc (size_t n, size_t sz)
{
  void * out;
  out = calloc (n, sz);
  if (!out) {
    fprintf(stderr,
	    "Out of memory trying to allocate/clear %zu items of %zu byte(s)\n",
	    n, sz);
    abort ();
  }
  return out;
}
void
xfree (void * p, size_t sz)
{
  free (p);
}
#endif

/* Spread the two 64-bit numbers into five nonzero values in the correct
 * range. */
void make_mrg_seed(uint64_t userseed1, uint64_t userseed2, uint_fast32_t* seed) {
  seed[0] = (uint32_t)(userseed1 & UINT32_C(0x3FFFFFFF)) + 1;
  seed[1] = (uint32_t)((userseed1 >> 30) & UINT32_C(0x3FFFFFFF)) + 1;
  seed[2] = (uint32_t)(userseed2 & UINT32_C(0x3FFFFFFF)) + 1;
  seed[3] = (uint32_t)((userseed2 >> 30) & UINT32_C(0x3FFFFFFF)) + 1;
  seed[4] = (uint32_t)((userseed2 >> 60) << 4) + (uint32_t)(userseed1 >> 60) + 1;
}
