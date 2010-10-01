/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <mpi.h>
#include <assert.h>
#include "common.h"

int rank, size;
#ifdef SIZE_MUST_BE_A_POWER_OF_TWO
int lgsize, size_minus_one;
#endif

void setup_globals() {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

#ifdef SIZE_MUST_BE_A_POWER_OF_TWO
  size_minus_one = size - 1;
  if (/* Check for power of 2 */ (size & (size - 1)) != 0) {
    fprintf(stderr, "Number of processes %d is not a power of two, yet SIZE_MUST_BE_A_POWER_OF_TWO is defined in main.cpp.\n", size);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  for (lgsize = 0; lgsize < size; ++lgsize) {
    if ((1 << lgsize) == size) break;
  }
  assert (lgsize < size);
#endif
}

void free_csr_graph(csr_graph* const g) {
  if (g->rowstarts != NULL) {free(g->rowstarts); g->rowstarts = NULL;}
  if (g->column != NULL) {free(g->column); g->column = NULL;}
}

/* These are in the graph generator. */
#if 0
void* xmalloc(size_t nbytes) {
  void* p = malloc(nbytes);
  if (!p) {
    fprintf(stderr, "malloc() failed for size %zu\n", nbytes);
    abort();
  }
  return p;
}

void* xcalloc(size_t n, size_t unit) {
  void* p = calloc(n, unit);
  if (!p) {
    fprintf(stderr, "calloc() failed for size %zu * %zu\n", n, unit);
    abort();
  }
  return p;
}
#endif

void* xrealloc(void* p, size_t nbytes) {
  p = realloc(p, nbytes);
  if (!p && nbytes != 0) {
    fprintf(stderr, "realloc() failed for size %zu\n", nbytes);
    abort();
  }
  return p;
}

void* xMPI_Alloc_mem(size_t nbytes) {
  void* p;
  MPI_Alloc_mem(nbytes, MPI_INFO_NULL, &p);
  if (!p) {
    fprintf(stderr, "MPI_Alloc_mem failed for size %zu\n", nbytes);
    abort();
  }
  return p;
}
