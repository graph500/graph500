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

static inline void* safe_malloc(size_t n) {
  void* p = malloc(n);
  if (!p) {
    fprintf(stderr, "Out of memory trying to allocate %zu byte(s)\n", n);
    abort();
  }
  return p;
}

static inline void* safe_calloc(size_t n, size_t k) {
  void* p = calloc(n, k);
  if (!p) {
    fprintf(stderr, "Out of memory trying to allocate %zu byte(s)\n", n);
    abort();
  }
  return p;
}

/* Spread the two 64-bit numbers into five nonzero values in the correct
 * range. */
static void make_mrg_seed(uint64_t userseed1, uint64_t userseed2, uint_fast32_t* seed) {
  seed[0] = (userseed1 & 0x3FFFFFFF) + 1;
  seed[1] = ((userseed1 >> 30) & 0x3FFFFFFF) + 1;
  seed[2] = (userseed2 & 0x3FFFFFFF) + 1;
  seed[3] = ((userseed2 >> 30) & 0x3FFFFFFF) + 1;
  seed[4] = ((userseed2 >> 60) << 4) + (userseed1 >> 60) + 1;
}

#ifdef GRAPH_GENERATOR_SEQ
void make_graph(int log_numverts, int64_t desired_nedges, uint64_t userseed1, uint64_t userseed2, double initiator[4], int64_t* nedges_ptr, int64_t** result_ptr) {
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
  generated_edge* edges = (generated_edge*)safe_calloc(nedges, sizeof(generated_edge)); /* multiplicity set to 0 for unused edges */
#else
  int64_t* edges = (int64_t*)safe_malloc(2 * nedges * sizeof(int64_t));
#endif

  generate_kronecker(0, 1, seed, log_numverts, M, initiator, edges);

  int64_t* vertex_perm = (int64_t*)safe_malloc(N * sizeof(int64_t));
  int64_t* result;
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  result = (int64_t*)safe_malloc(2 * nedges * sizeof(int64_t));
#else
  result = edges;
#endif
  *result_ptr = result;

  mrg_state state;
  mrg_seed(&state, seed[0], seed[1], seed[2], seed[3], seed[4]);
  rand_sort_shared(&state, N, vertex_perm);
  int64_t i;
  /* Apply vertex permutation to graph, optionally copying into user's result
   * array. */
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  for (i = 0; i < nedges; ++i) {
    if (edges[i].multiplicity != 0) {
      result[i * 2] = vertex_perm[edges[i].src];
      result[i * 2 + 1] = vertex_perm[edges[i].tgt];
    } else {
      result[i * 2] = result[i * 2 + 1] = (int64_t)(-1);
    }
  }
  free(edges);
#else
  for (i = 0; i < 2 * nedges; ++i) {
    if (edges[i] != (int64_t)(-1)) {
      edges[i] = vertex_perm[edges[i]];
    }
  }
#endif

  free(vertex_perm);
}
#endif /* GRAPH_GENERATOR_SEQ */

#ifdef __MTA__
void make_graph(int log_numverts, int64_t desired_nedges, uint64_t userseed1, uint64_t userseed2, double initiator[4], int64_t* nedges_ptr, int64_t** result_ptr) {
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
  generated_edge* edges = (generated_edge*)safe_calloc(nedges, sizeof(generated_edge)); /* multiplicity set to 0 for unused edges */
#else
  int64_t* edges = (int64_t*)safe_malloc(2 * nedges * sizeof(int64_t));
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

  int64_t* vertex_perm = (int64_t*)safe_malloc(N * sizeof(int64_t));
  int64_t* result;

#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  result = (int64_t*)safe_malloc(2 * nedges * sizeof(int64_t));
#else
  result = edges;
#endif
  *result_ptr = result;

  mrg_state state;
  mrg_seed(&state, seed[0], seed[1], seed[2], seed[3], seed[4]);
  rand_sort_shared(&state, N, vertex_perm);
  int64_t i;
  /* Apply vertex permutation to graph, optionally copying into user's result
   * array. */
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
#pragma mta assert parallel
#pragma mta block schedule
  for (i = 0; i < nedges; ++i) {
    if (edges[i].multiplicity != 0) {
      result[i * 2] = vertex_perm[edges[i].src];
      result[i * 2 + 1] = vertex_perm[edges[i].tgt];
    } else {
      result[i * 2] = result[i * 2 + 1] = (int64_t)(-1);
    }
  }
  free(edges);
#else
#pragma mta assert parallel
#pragma mta block schedule
  for (i = 0; i < 2 * nedges; ++i) {
    if (edges[i] != (int64_t)(-1)) {
      edges[i] = vertex_perm[edges[i]];
    }
  }
#endif

  free(vertex_perm);
}
#endif /* __MTA__ */

#ifdef GRAPH_GENERATOR_OMP
void make_graph(int log_numverts, int64_t desired_nedges, uint64_t userseed1, uint64_t userseed2, double initiator[4], int64_t* nedges_ptr, int64_t** result_ptr) {
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
  generated_edge* edges = (generated_edge*)safe_calloc(nedges, sizeof(generated_edge)); /* multiplicity set to 0 for unused edges */
#else
  int64_t* edges = (int64_t*)safe_malloc(2 * nedges * sizeof(int64_t));
#endif

#pragma omp parallel
  {
    int rank = omp_get_thread_num(), size = omp_get_num_threads();
    generate_kronecker(rank, size, seed, log_numverts, M, initiator, edges);
  }

  int64_t* vertex_perm = (int64_t*)safe_malloc(N * sizeof(int64_t));
  int64_t* result;

#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  result = (int64_t*)safe_malloc(2 * nedges * sizeof(int64_t));
#else
  result = edges;
#endif
  *result_ptr = result;

  mrg_state state;
  mrg_seed(&state, seed[0], seed[1], seed[2], seed[3], seed[4]);
  rand_sort_shared(&state, N, vertex_perm);
  int64_t i;
  /* Apply vertex permutation to graph, optionally copying into user's result
   * array. */
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
#pragma omp parallel for
  for (i = 0; i < nedges; ++i) {
    if (edges[i].multiplicity != 0) {
      result[i * 2] = vertex_perm[edges[i].src];
      result[i * 2 + 1] = vertex_perm[edges[i].tgt];
    } else {
      result[i * 2] = result[i * 2 + 1] = (int64_t)(-1);
    }
  }
  free(edges);
#else
#pragma omp parallel for
  for (i = 0; i < 2 * nedges; ++i) {
    if (edges[i] != (int64_t)(-1)) {
      edges[i] = vertex_perm[edges[i]];
    }
  }
#endif

  free(vertex_perm);
}
#endif /* GRAPH_GENERATOR_OMP */

#ifdef GRAPH_GENERATOR_MPI
static inline int compare_int64(const void* a, const void* b) {
  int64_t aa = *(const int64_t*)a;
  int64_t bb = *(const int64_t*)b;
  return (aa < bb) ? -1 : (aa == bb) ? 0 : 1;
}

void make_graph(int log_numverts, int64_t desired_nedges, uint64_t userseed1, uint64_t userseed2, double initiator[4], int64_t* nedges_ptr, int64_t** result_ptr) {
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
  *nedges_ptr = nedges;
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  generated_edge* local_edges = (generated_edge*)safe_calloc(nedges, sizeof(generated_edge)); /* multiplicity set to 0 for unused edges */
#else
  int64_t* local_edges = (int64_t*)safe_malloc(2 * nedges * sizeof(int64_t));
#endif

  generate_kronecker(rank, size, seed, log_numverts, M, initiator, local_edges);

  int64_t* local_vertex_perm = NULL;

  mrg_state state;
  mrg_seed(&state, seed[0], seed[1], seed[2], seed[3], seed[4]);
  int64_t perm_local_size;
  rand_sort_mpi(MPI_COMM_WORLD, &state, N, &perm_local_size, &local_vertex_perm);

  /* Apply vertex permutation to graph, optionally copying into user's result
   * array. */

  /* Get the size of the permutation fragment on each node. */
  int64_t* perm_sizes = (int64_t*)safe_malloc(size * sizeof(int64_t));
  MPI_Allgather(&perm_local_size, 1, INT64_T_MPI_TYPE, perm_sizes, 1, INT64_T_MPI_TYPE, MPI_COMM_WORLD);

  /* Make a table consisting of the owner of each multiple of minpermsize (the
   * minimum permutation fragment size across all nodes).  Using this table,
   * only two locations need to be considered to find the owner of any given
   * vertex in the distributed permutation, avoiding a binary search.  Also,
   * sanity check the permutation to make sure its distribution is reasonable
   * (PRNG problems can cause bad distributions). */
  int minpermsize = perm_local_size, maxpermsize = perm_local_size;
  int64_t i;
  for (i = 0; i < size; ++i) {
    if (perm_sizes[i] < minpermsize) minpermsize = perm_sizes[i];
    if (perm_sizes[i] > maxpermsize) maxpermsize = perm_sizes[i];
  }

  if (minpermsize == 0) {
    fprintf(stderr, "One node has permutation size 0\n");
    abort();
  }
  if (maxpermsize > 2 * minpermsize && rank == 0) {
    fprintf(stderr, "Minimum permutation bucket size %d is too small compared to maximum size %d\n", minpermsize, maxpermsize);
    abort();
  }

  /* Prefix sum permutation sizes (exclusive) for owner lookups. */
  int* perm_displs = (int*)safe_malloc((size + 1) * sizeof(int));
  perm_displs[0] = 0;
  for (i = 1; i < size + 1; ++i) perm_displs[i] = perm_displs[i - 1] + perm_sizes[i - 1];
  int64_t my_perm_displ = perm_displs[rank];
  free(perm_sizes); perm_sizes = NULL;

  /* Copy the edge endpoints into the result array if necessary. */
  int64_t* result;
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  result = (int64_t*)safe_malloc(2 * nedges * sizeof(int64_t));
  for (i = 0; i < nedges; ++i) {
    if (local_edges[i].multiplicity != 0) {
      result[i * 2] = local_edges[i].src;
      result[i * 2 + 1] = local_edges[i].tgt;
    } else {
      result[i * 2] = result[i * 2 + 1] = (int64_t)(-1);
    }
  }
  free(local_edges); local_edges = NULL;
  *result_ptr = result;
#else
  result = local_edges;
  *result_ptr = result;
  local_edges = NULL; /* Freed by caller */
#endif

  /* Create a sorted list, without duplicates, of all vertices mentioned in the
   * local edge list. */
  int64_t* vertices_used = (int64_t*)safe_malloc(2 * nedges * sizeof(int64_t));
  memcpy(vertices_used, result, 2 * nedges * sizeof(int64_t));

  /* FIXME: use better sorting function */
  qsort(vertices_used, 2 * nedges, sizeof(int64_t), &compare_int64);
  int64_t n_unique_vertices;
  {
    /* Remove duplicates in-place. */
    int64_t i_in, i_out;
    for (i_in = 0, i_out = 0; i_in < 2 * nedges; ++i_in) {
      if (vertices_used[i_in] != -1 &&
          (i_in == 0 || vertices_used[i_in] != vertices_used[i_in - 1])) {
        vertices_used[i_out++] = vertices_used[i_in];
      } else { /* Skip this input element. */ }
    }
    /* Shrink buffer. */
    vertices_used = realloc(vertices_used, i_out * sizeof(int64_t));
    if (!vertices_used) {
      fprintf(stderr, "Shrinking buffer failed; should not happen.\n");
      abort();
    }
    n_unique_vertices = i_out;
  }

  int64_t* lookup_results = (int64_t*)safe_malloc(n_unique_vertices * sizeof(int64_t));

  /* Send the set of unique vertices owned by a given processor to that
   * processor, map them to the correct permutation entries there, and send the
   * replies back. */
  {
    int64_t* msg_locs = (int64_t*)safe_malloc((size + 1) * sizeof(int64_t));
    msg_locs[0] = 0;
    int dest = 1;
    for (i = 0; i < n_unique_vertices; ++i) {
      while (dest <= size && vertices_used[i] >= perm_displs[dest]) {
        msg_locs[dest] = i;
        ++dest;
      }
    }
    for (; dest <= size; ++dest) {
      msg_locs[dest] = n_unique_vertices;
    }
    free(perm_displs); perm_displs = NULL;

    int* will_send_request = (int*)safe_malloc(size * sizeof(int));
    for (dest = 0; dest < size; ++dest) {
      will_send_request[dest] = (msg_locs[dest] != msg_locs[dest + 1]) && (dest != rank);
    }
    int num_reqs_to_receive;
    int* ones = (int*)safe_malloc(size * sizeof(int)); /* Data sizes */
    for (i = 0; i < size; ++i) ones[i] = 1;
    MPI_Reduce_scatter(will_send_request, &num_reqs_to_receive, ones, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    free(will_send_request); will_send_request = NULL;
    free(ones); ones = NULL;

    /* Do our own lookups locally. */
    for (i = msg_locs[rank]; i < msg_locs[rank + 1]; ++i) {
      lookup_results[i] = local_vertex_perm[vertices_used[i] - my_perm_displ];
    }

    /* Have only one request of each type outstanding at a time to save memory.
     * Requests use ssend as a cheap way to do flow control; replies use rsend
     * to allow eager protocols to be used more often. */
    int destoffset = 1;
    MPI_Request sendreq, recvreq, response_recvreq;
    int send_pending = 0, recv_pending = 0;
    int64_t* incoming_lookups = (int64_t*)safe_malloc(perm_local_size * sizeof(int64_t));
    int64_t* incoming_responses = (int64_t*)safe_malloc(perm_local_size * sizeof(int64_t));
    while (1) {
      if (!send_pending && destoffset < size) {
        int dest = (rank + destoffset) % size;
        if (msg_locs[dest] == msg_locs[dest + 1]) {
          ++destoffset;
        } else {
          /* Recv must be first to allow rsend on receiver. */
          MPI_Irecv(&lookup_results[msg_locs[dest]], msg_locs[dest + 1] - msg_locs[dest],
                    INT64_T_MPI_TYPE, dest, 1, MPI_COMM_WORLD, &response_recvreq);
          MPI_Issend(&vertices_used[msg_locs[dest]], msg_locs[dest + 1] - msg_locs[dest],
                     INT64_T_MPI_TYPE, dest, 0, MPI_COMM_WORLD, &sendreq);
          send_pending = 1;
        }
      }
      if (!recv_pending && num_reqs_to_receive > 0) {
        MPI_Irecv(incoming_lookups, perm_local_size, INT64_T_MPI_TYPE,
                  MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recvreq);
        --num_reqs_to_receive;
        recv_pending = 1;
      }
      if (send_pending) {
        int flag;
        MPI_Test(&sendreq, &flag, MPI_STATUS_IGNORE);
        if (flag) {
          send_pending = 0;
          MPI_Wait(&response_recvreq, MPI_STATUS_IGNORE);
          ++destoffset;
        }
      }
      if (recv_pending) {
        int flag;
        MPI_Status st;
        MPI_Test(&recvreq, &flag, &st);
        if (flag) {
          int nlookups;
          MPI_Get_count(&st, INT64_T_MPI_TYPE, &nlookups);
          for (i = 0; i < nlookups; ++i) {
            assert (incoming_lookups[i] >= my_perm_displ &&
                    incoming_lookups[i] < my_perm_displ + perm_local_size);
            incoming_responses[i] = local_vertex_perm[incoming_lookups[i] - my_perm_displ];
          }
          MPI_Rsend(incoming_responses, nlookups, INT64_T_MPI_TYPE,
                    st.MPI_SOURCE, 1, MPI_COMM_WORLD);
          recv_pending = 0;
        }
      }
      if (!send_pending && !recv_pending &&
          destoffset >= size && num_reqs_to_receive == 0) {
        break;
      }
    }
    free(msg_locs); msg_locs = NULL;
    free(incoming_lookups); incoming_lookups = NULL;
    free(incoming_responses); incoming_responses = NULL;
  }

  /* Do the edge endpoint permutation. */
  for (i = 0; i < 2 * nedges; ++i) {
    int64_t orig_vertex = result[i];
    if (orig_vertex != -1) {
      /* Binary search in zip(vertices_used, lookup_results). */
      /* FIXME: See if a hash table or such would be better. */
      int64_t low = 0 /* Inclusive */,
              high = n_unique_vertices /* Exclusive */;
      while (high > low + 1) {
        int64_t mid = low + (high - low) / 2;
        if (orig_vertex == vertices_used[mid]) {
          low = mid; /* Found */
          high = low + 1;
          break;
        } else if (orig_vertex > vertices_used[mid]) {
          low = mid + 1;
          if (low == high) break; /* Not found */
        } else {
          high = mid;
        }
      }
      assert (high == low + 1);
      assert (vertices_used[low] == orig_vertex);
      result[i] = lookup_results[low];
    }
  }

  free(vertices_used); vertices_used = NULL;
  free(lookup_results); lookup_results = NULL;
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
  mrg_seed(&st, seed[0], seed[1], seed[2], seed[3], seed[4]);

  mrg_skip(&st, 2, 0, 2 * position); /* Each double takes two PRNG outputs */

  for (i = 0; i < nvalues; ++i) {
    result[i] = mrg_get_double_orig(&st);
  }
}
