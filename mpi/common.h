/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <stddef.h>

#define INT64_T_MPI_TYPE MPI_LONG_LONG

#define SIZE_MUST_BE_A_POWER_OF_TWO

extern int rank, size;
#ifdef SIZE_MUST_BE_A_POWER_OF_TWO
extern int lgsize, size_minus_one;
#endif

/* Distribute edges by their endpoints (make two directed copies of each input
 * undirected edge); distribution is 1-d and cyclic. */
#ifdef SIZE_MUST_BE_A_POWER_OF_TWO
#define MOD_SIZE(v) ((v) & size_minus_one)
#define DIV_SIZE(v) ((v) >> lgsize)
#else
#define MOD_SIZE(v) ((v) % size)
#define DIV_SIZE(v) ((v) / size)
#endif
#define VERTEX_OWNER(v) ((int)(MOD_SIZE(v)))
#define VERTEX_LOCAL(v) ((size_t)(DIV_SIZE(v)))
#define VERTEX_TO_GLOBAL(i) ((int64_t)((i) * size + rank))

typedef struct csr_graph {
  size_t nlocalverts;
  size_t nlocaledges;
  int64_t nglobalverts;
  size_t *rowstarts;
  int64_t *column;
} csr_graph;

#ifdef __cplusplus
extern "C" {
#endif

void setup_globals(void); /* In utils.c */
void free_csr_graph(csr_graph* const g); /* In utils.c */
void* xMPI_Alloc_mem(size_t nbytes); /* In utils.c */
void* xmalloc(size_t nbytes); /* In utils.c */
void* xcalloc(size_t n, size_t unit); /* In utils.c */
void* xrealloc(void* p, size_t nbytes); /* In utils.c */

void convert_graph_to_csr(const int64_t nedges, const int64_t* const edges, csr_graph* const g); /* In convert_to_csr.c */
void find_bfs_roots(int *num_bfs_roots, const csr_graph* const g, const uint64_t seed1, const uint64_t seed2, int64_t* const bfs_roots); /* In find_roots.c */
int validate_bfs_result(const csr_graph* const g, const int64_t root, const int64_t* const pred, const int64_t nvisited); /* In validate.c */

void run_mpi_bfs(const csr_graph* const g, int64_t root, int64_t* pred, int64_t* nvisited); /* Provided by user */

#ifdef __cplusplus
}
#endif

#endif /* COMMON_H */
