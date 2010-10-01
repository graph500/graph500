/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#ifdef __MTA__
#include <sys/mta_task.h>
#endif
#ifdef GRAPH_GENERATOR_MPI
#include <mpi.h>
#endif
#ifdef GRAPH_GENERATOR_OMP
#include <omp.h>
#endif

/* Simplified interface to build graphs with scrambled vertices. */

#include "graph_generator.h"
#include "permutation_gen.h"
#include "apply_permutation_mpi.h"
#include "scramble_edges.h"
#include "utils.h"

#ifdef GRAPH_GENERATOR_SEQ
void make_graph(int log_numverts, int64_t desired_nedges, uint64_t userseed1, uint64_t userseed2, const double initiator[4], int64_t* nedges_ptr, int64_t** result_ptr) {
  int64_t N, M;

  N = (int64_t)pow(GRAPHGEN_INITIATOR_SIZE, log_numverts);
  M = desired_nedges;

  /* Spread the two 64-bit numbers into five nonzero values in the correct
   * range. */
  uint_fast32_t seed[5];
  make_mrg_seed(userseed1, userseed2, seed);

  int64_t nedges = compute_edge_array_size(0, 1, M);
  *nedges_ptr = nedges;
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  generated_edge* edges = (generated_edge*)xcalloc(nedges, sizeof(generated_edge)); /* multiplicity set to 0 for unused edges */
#else
  int64_t* edges = (int64_t*)xmalloc(2 * nedges * sizeof(int64_t));
#endif

  generate_kronecker(0, 1, seed, log_numverts, M, initiator, edges);

  int64_t* vertex_perm = (int64_t*)xmalloc(N * sizeof(int64_t));
  int64_t* result;
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  result = (int64_t*)xmalloc(2 * nedges * sizeof(int64_t));
#else
  result = edges;
#endif
  *result_ptr = result;

  mrg_state state;
  mrg_seed(&state, seed);
  rand_sort_shared(&state, N, vertex_perm);
  int64_t i;
  /* Apply vertex permutation to graph, optionally copying into user's result
   * array. */
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  for (i = 0; i < nedges; ++i) {
    if (edges[i].multiplicity != 0) {
      int64_t v1 = vertex_perm[edges[i].src];
      int64_t v2 = vertex_perm[edges[i].tgt];
      /* Sort these since otherwise the directions of the permuted edges would
       * give away the unscrambled vertex order. */
      result[i * 2] = (v1 < v2) ? v1 : v2;
      result[i * 2 + 1] = (v1 < v2) ? v2 : v1;
    } else {
      result[i * 2] = result[i * 2 + 1] = (int64_t)(-1);
    }
  }
  free(edges);
#else
  for (i = 0; i < 2 * nedges; i += 2) {
    if (edges[i] != (int64_t)(-1)) {
      int64_t v1 = vertex_perm[edges[i]];
      int64_t v2 = vertex_perm[edges[i + 1]];
      /* Sort these since otherwise the directions of the permuted edges would
       * give away the unscrambled vertex order. */
      edges[i] = (v1 < v2) ? v1 : v2;
      edges[i + 1] = (v1 < v2) ? v2 : v1;
    }
  }
#endif

  free(vertex_perm);

  /* Randomly mix up the order of the edges. */
  scramble_edges_shared(userseed1, userseed2, nedges, edges);
}
#endif /* GRAPH_GENERATOR_SEQ */

#ifdef __MTA__
void make_graph(int log_numverts, int64_t desired_nedges, uint64_t userseed1, uint64_t userseed2, const double initiator[4], int64_t* nedges_ptr, int64_t** result_ptr) {
  int64_t N, M;

  N = (int64_t)pow(GRAPHGEN_INITIATOR_SIZE, log_numverts);
  M = (int64_t)desired_nedges;

  /* Spread the two 64-bit numbers into five nonzero values in the correct
   * range. */
  uint_fast32_t seed[5];
  make_mrg_seed(userseed1, userseed2, seed);

  int64_t nedges = compute_edge_array_size(0, 1, M);
  *nedges_ptr = nedges;
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  generated_edge* edges = (generated_edge*)xcalloc(nedges, sizeof(generated_edge)); /* multiplicity set to 0 for unused edges */
#else
  int64_t* edges = (int64_t*)xmalloc(2 * nedges * sizeof(int64_t));
#endif

  int rank, size;
  /* The "for all streams" here is in compiler versions >= 6.4 */
#pragma mta use 100 streams
#pragma mta for all streams rank of size
  {
    double my_initiator[GRAPHGEN_INITIATOR_SIZE * GRAPHGEN_INITIATOR_SIZE]; /* Local copy */
    int i;
    for (i = 0; i < GRAPHGEN_INITIATOR_SIZE * GRAPHGEN_INITIATOR_SIZE; ++i) {
      my_initiator[i] = initiator[i];
    }
    generate_kronecker(rank, size, seed, log_numverts, M, my_initiator, edges);
  }

  int64_t* vertex_perm = (int64_t*)xmalloc(N * sizeof(int64_t));
  int64_t* result;

#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  result = (int64_t*)xmalloc(2 * nedges * sizeof(int64_t));
#else
  result = edges;
#endif
  *result_ptr = result;

  mrg_state state;
  mrg_seed(&state, seed);
  rand_sort_shared(&state, N, vertex_perm);
  int64_t i;
  /* Apply vertex permutation to graph, optionally copying into user's result
   * array. */
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
#pragma mta assert parallel
#pragma mta block schedule
  for (i = 0; i < nedges; ++i) {
    if (edges[i].multiplicity != 0) {
      int64_t v1 = vertex_perm[edges[i].src];
      int64_t v2 = vertex_perm[edges[i].tgt];
      /* Sort these since otherwise the directions of the permuted edges would
       * give away the unscrambled vertex order. */
      result[i * 2] = MTA_INT_MIN(v1, v2);
      result[i * 2 + 1] = MTA_INT_MAX(v1, v2);
    } else {
      result[i * 2] = result[i * 2 + 1] = (int64_t)(-1);
    }
  }
  free(edges);
#else
#pragma mta assert parallel
#pragma mta block schedule
  for (i = 0; i < 2 * nedges; i += 2) {
    if (edges[i] != (int64_t)(-1)) {
      int64_t v1 = vertex_perm[edges[i]];
      int64_t v2 = vertex_perm[edges[i + 1]];
      /* Sort these since otherwise the directions of the permuted edges would
       * give away the unscrambled vertex order. */
      edges[i] = MTA_INT_MIN(v1, v2);
      edges[i + 1] = MTA_INT_MAX(v1, v2);
    }
  }
#endif

  free(vertex_perm);

  /* Randomly mix up the order of the edges. */
  scramble_edges_shared(userseed1, userseed2, nedges, edges);
}
#endif /* __MTA__ */

#ifdef GRAPH_GENERATOR_OMP
void make_graph(int log_numverts, int64_t desired_nedges, uint64_t userseed1, uint64_t userseed2, const double initiator[4], int64_t* nedges_ptr, int64_t** result_ptr) {
  int64_t N, M;

  N = (int64_t)pow(GRAPHGEN_INITIATOR_SIZE, log_numverts);
  M = desired_nedges;

  /* Spread the two 64-bit numbers into five nonzero values in the correct
   * range. */
  uint_fast32_t seed[5];
  make_mrg_seed(userseed1, userseed2, seed);

  int64_t nedges = compute_edge_array_size(0, 1, M);
  *nedges_ptr = nedges;
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  generated_edge* edges = (generated_edge*)xcalloc(nedges, sizeof(generated_edge)); /* multiplicity set to 0 for unused edges */
#else
  int64_t* edges = (int64_t*)xmalloc(2 * nedges * sizeof(int64_t));
#endif

#pragma omp parallel
  {
    int rank = omp_get_thread_num(), size = omp_get_num_threads();
    generate_kronecker(rank, size, seed, log_numverts, M, initiator, edges);
  }

  int64_t* vertex_perm = (int64_t*)xmalloc(N * sizeof(int64_t));
  int64_t* result;

#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  result = (int64_t*)xmalloc(2 * nedges * sizeof(int64_t));
#else
  result = edges;
#endif
  *result_ptr = result;

  mrg_state state;
  mrg_seed(&state, seed);
  rand_sort_shared(&state, N, vertex_perm);
  int64_t i;
  /* Apply vertex permutation to graph, optionally copying into user's result
   * array. */
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
#pragma omp parallel for
  for (i = 0; i < nedges; ++i) {
    if (edges[i].multiplicity != 0) {
      int64_t v1 = vertex_perm[edges[i].src];
      int64_t v2 = vertex_perm[edges[i].tgt];
      /* Sort these since otherwise the directions of the permuted edges would
       * give away the unscrambled vertex order. */
      result[i * 2] = (v1 < v2) ? v1 : v2;
      result[i * 2 + 1] = (v1 < v2) ? v2 : v1;
    } else {
      result[i * 2] = result[i * 2 + 1] = (int64_t)(-1);
    }
  }
  free(edges);
#else
#pragma omp parallel for
  for (i = 0; i < 2 * nedges; i += 2) {
    if (edges[i] != (int64_t)(-1)) {
      int64_t v1 = vertex_perm[edges[i]];
      int64_t v2 = vertex_perm[edges[i + 1]];
      /* Sort these since otherwise the directions of the permuted edges would
       * give away the unscrambled vertex order. */
      edges[i] = (v1 < v2) ? v1 : v2;
      edges[i + 1] = (v1 < v2) ? v2 : v1;
    }
  }
#endif

  free(vertex_perm);

  /* Randomly mix up the order of the edges. */
  scramble_edges_shared(userseed1, userseed2, nedges, edges);
}
#endif /* GRAPH_GENERATOR_OMP */

#ifdef GRAPH_GENERATOR_MPI
void make_graph(int log_numverts, int64_t desired_nedges, uint64_t userseed1, uint64_t userseed2, const double initiator[4], int64_t* nedges_ptr, int64_t** result_ptr) {
  int64_t N, M;
  int rank, size;

  N = (int64_t)pow(GRAPHGEN_INITIATOR_SIZE, log_numverts);
  M = desired_nedges;

  /* Spread the two 64-bit numbers into five nonzero values in the correct
   * range. */
  uint_fast32_t seed[5];
  make_mrg_seed(userseed1, userseed2, seed);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int64_t nedges = compute_edge_array_size(rank, size, M);
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  generated_edge* local_edges = (generated_edge*)xmalloc(nedges * sizeof(generated_edge));
#else
  int64_t* local_edges = (int64_t*)xmalloc(2 * nedges * sizeof(int64_t));
#endif

  double start = MPI_Wtime();
  generate_kronecker(rank, size, seed, log_numverts, M, initiator, local_edges);
  double gen_time = MPI_Wtime() - start;

  int64_t* local_vertex_perm = NULL;

  mrg_state state;
  mrg_seed(&state, seed);
  start = MPI_Wtime();
  int64_t perm_local_size;
  rand_sort_mpi(MPI_COMM_WORLD, &state, N, &perm_local_size, &local_vertex_perm);
  double perm_gen_time = MPI_Wtime() - start;

  /* Copy the edge endpoints into the result array if necessary. */
  int64_t* result;
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  result = (int64_t*)xmalloc(2 * nedges * sizeof(int64_t));
  for (i = 0; i < nedges; ++i) {
    if (local_edges[i].multiplicity != 0) {
      result[i * 2] = local_edges[i].src;
      result[i * 2 + 1] = local_edges[i].tgt;
    } else {
      result[i * 2] = result[i * 2 + 1] = (int64_t)(-1);
    }
  }
  free(local_edges); local_edges = NULL;
#else
  result = local_edges;
  *result_ptr = result;
  local_edges = NULL; /* Freed by caller */
#endif

  /* Apply vertex permutation to graph. */
  start = MPI_Wtime();
  apply_permutation_mpi(MPI_COMM_WORLD, perm_local_size, local_vertex_perm, N, nedges, result);
  double perm_apply_time = MPI_Wtime() - start;

  free(local_vertex_perm); local_vertex_perm = NULL;

  /* Randomly mix up the order of the edges. */
  start = MPI_Wtime();
  int64_t* new_result;
  int64_t nedges_out;
  scramble_edges_mpi(MPI_COMM_WORLD, userseed1, userseed2, nedges, result, &nedges_out, &new_result);
  double edge_scramble_time = MPI_Wtime() - start;

  free(result); result = NULL;

  *result_ptr = new_result;
  *nedges_ptr = nedges_out;

  if (rank == 0) {
    fprintf(stdout, "unpermuted_graph_generation:    %f s\n", gen_time);
    fprintf(stdout, "vertex_permutation_generation:  %f s\n", perm_gen_time);
    fprintf(stdout, "vertex_permutation_application: %f s\n", perm_apply_time);
    fprintf(stdout, "edge_scrambling:                %f s\n", edge_scramble_time);
  }
}
#endif

/* PRNG interface for implementations; takes seed in same format as given by
 * users, and creates a vector of doubles in a reproducible (and
 * random-access) way. */
void make_random_numbers(
       /* in */ int64_t nvalues    /* Number of values to generate */,
       /* in */ uint64_t userseed1 /* Arbitrary 64-bit seed value */,
       /* in */ uint64_t userseed2 /* Arbitrary 64-bit seed value */,
       /* in */ int64_t position   /* Start index in random number stream */,
       /* out */ double* result    /* Returned array of values */
) {
  int64_t i;
  uint_fast32_t seed[5];
  make_mrg_seed(userseed1, userseed2, seed);

  mrg_state st;
  mrg_seed(&st, seed);

  mrg_skip(&st, 2, 0, 2 * position); /* Each double takes two PRNG outputs */

  for (i = 0; i < nvalues; ++i) {
    result[i] = mrg_get_double_orig(&st);
  }
}
