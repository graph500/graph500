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
#ifdef __MTA__
#include <sys/mta_task.h>
#endif
#ifdef GRAPH_GENERATOR_MPI
#include <mpi.h>
#endif

/* Simplified interface for build graphs with scrambled vertices. */

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

#ifdef GRAPH_GENERATOR_SEQ
void make_graph(int log_numverts, uint64_t userseed1, uint64_t userseed2, count_type* nedges, count_type** result) {
  count_type N, M;

  double avg_degree = 16.;
#ifdef GRAPHGEN_UNDIRECTED
  avg_degree /= 2; /* Adjust for edge doubling */
#endif
  N = count_pow(2, log_numverts);
  M = (count_type)(N * avg_degree);

  /* Spread the two 64-bit numbers into five nonzero values in the correct
   * range. */
  uint_fast32_t seed[5];
  seed[0] = (userseed1 & 0x3FFFFFFF) + 1;
  seed[1] = ((userseed1 >> 30) & 0x3FFFFFFF) + 1;
  seed[2] = (userseed2 & 0x3FFFFFFF) + 1;
  seed[3] = ((userseed2 >> 30) & 0x3FFFFFFF) + 1;
  seed[4] = ((userseed2 >> 60) << 4) + (userseed1 >> 60) + 1;

  *nedges = compute_edge_array_size(0, 1, M);
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  generated_edge* edges = (generated_edge*)safe_calloc(*nedges, sizeof(generated_edge)); /* multiplicity set to 0 for unused edges */
#else
  count_type* edges = (count_type*)safe_malloc(2 * *nedges * sizeof(count_type));
#endif

  double initiator[2 * 2] = {.57, .19, .19, .05};
  generate_kronecker(0, 1, seed, log_numverts, M, initiator, edges);

  count_type* vertex_perm = (count_type*)safe_malloc(N * sizeof(count_type));
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  *result = (count_type*)safe_malloc(2 * *nedges * sizeof(count_type));
#else
  *result = edges;
#endif

  mrg_state state;
  mrg_seed(&state, seed[0], seed[1], seed[2], seed[3], seed[4]);
  rand_sort_shared(&state, N, vertex_perm);
  count_type i;
  /* Apply vertex permutation to graph, optionally copying into user's result
   * array. */
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  for (i = 0; i < *nedges; ++i) {
    if (edges[i].multiplicity != 0) {
      (*result)[i * 2] = vertex_perm[edges[i].src];
      (*result)[i * 2 + 1] = vertex_perm[edges[i].tgt];
    } else {
      (*result)[i * 2] = (*result)[i * 2 + 1] = (count_type)(-1);
    }
  }
  free(edges);
#else
  for (i = 0; i < 2 * *nedges; ++i) {
    if (edges[i] != (count_type)(-1)) {
      edges[i] = vertex_perm[edges[i]];
    }
  }
#endif

  free(vertex_perm);
}
#endif /* GRAPH_GENERATOR_SEQ */

#ifdef __MTA__
void make_graph(int log_numverts, uint64_t userseed1, uint64_t userseed2, count_type* nedges, count_type** result) {
  count_type N, M;

  double avg_degree = 16.;
#ifdef GRAPHGEN_UNDIRECTED
  avg_degree /= 2; /* Adjust for edge doubling */
#endif
  N = count_pow(2, log_numverts);
  M = (count_type)(N * avg_degree);

  /* Spread the two 64-bit numbers into five nonzero values in the correct
   * range. */
  uint_fast32_t seed[5];
  seed[0] = (userseed1 & 0x3FFFFFFF) + 1;
  seed[1] = ((userseed1 >> 30) & 0x3FFFFFFF) + 1;
  seed[2] = (userseed2 & 0x3FFFFFFF) + 1;
  seed[3] = ((userseed2 >> 30) & 0x3FFFFFFF) + 1;
  seed[4] = ((userseed2 >> 60) << 4) + (userseed1 >> 60) + 1;

  *nedges = compute_edge_array_size(0, 1, M);
  generated_edge* local_edges = (generated_edge*)safe_calloc(*nedges, sizeof(generated_edge)); /* multiplicity set to 0 for unused edges */

  int rank, size;
  /* The "for all streams" here is in compiler versions >= 6.4 */
#pragma mta use 100 streams
#pragma mta for all streams rank of size
  {
    double my_initiator[2 * 2] = {.57, .19, .19, .05}; /* Local copy */
    generate_kronecker(rank, size, seed, log_numverts, M, my_initiator, local_edges);
  }

  count_type* vertex_perm = (count_type*)safe_malloc(N * sizeof(count_type));

#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  *result = (count_type*)safe_malloc(2 * *nedges * sizeof(count_type));
#else
  *result = local_edges;
#endif

  mrg_state state;
  mrg_seed(&state, seed[0], seed[1], seed[2], seed[3], seed[4]);
  rand_sort_shared(&state, N, vertex_perm);
  count_type i;
  /* Apply vertex permutation to graph, optionally copying into user's result
   * array. */
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
#pragma mta assert parallel
#pragma mta block schedule
  for (i = 0; i < *nedges; ++i) {
    if (local_edges[i].multiplicity != 0) {
      (*result)[i * 2] = vertex_perm[local_edges[i].src];
      (*result)[i * 2 + 1] = vertex_perm[local_edges[i].tgt];
    } else {
      (*result)[i * 2] = (*result)[i * 2 + 1] = (count_type)(-1);
    }
  }
  free(local_edges);
#else
#pragma mta assert parallel
#pragma mta block schedule
  for (i = 0; i < 2 * *nedges; ++i) {
    if (local_edges[i] != (count_type)(-1)) {
      local_edges[i] = vertex_perm[local_edges[i]];
    }
  }
#endif

  free(vertex_perm);
}
#endif /* __MTA__ */

#ifdef GRAPH_GENERATOR_MPI
void make_graph(int log_numverts, uint64_t userseed1, uint64_t userseed2, count_type* nedges, count_type** result) {
  count_type N, M;
  int rank, size;

  double avg_degree = 16.;
#ifdef GRAPHGEN_UNDIRECTED
  avg_degree /= 2; /* Adjust for edge doubling */
#endif
  N = count_pow(2, log_numverts);
  M = (count_type)(N * avg_degree);

  /* Spread the two 64-bit numbers into five nonzero values in the correct
   * range. */
  uint_fast32_t seed[5];
  seed[0] = (userseed1 & 0x3FFFFFFF) + 1;
  seed[1] = ((userseed1 >> 30) & 0x3FFFFFFF) + 1;
  seed[2] = (userseed2 & 0x3FFFFFFF) + 1;
  seed[3] = ((userseed2 >> 30) & 0x3FFFFFFF) + 1;
  seed[4] = ((userseed2 >> 60) << 4) + (userseed1 >> 60) + 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  *nedges = compute_edge_array_size(rank, size, M);
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  generated_edge* local_edges = (generated_edge*)safe_calloc(*nedges, sizeof(generated_edge)); /* multiplicity set to 0 for unused edges */
#else
  count_type* local_edges = (count_type*)safe_malloc(2 * *nedges * sizeof(count_type));
#endif

  double initiator[2 * 2] = {.57, .19, .19, .05};
  generate_kronecker(rank, size, seed, log_numverts, M, initiator, local_edges);

  count_type* local_vertex_perm = NULL;

  mrg_state state;
  mrg_seed(&state, seed[0], seed[1], seed[2], seed[3], seed[4]);
  count_type perm_local_size;
  rand_sort_mpi(MPI_COMM_WORLD, &state, N, &perm_local_size, &local_vertex_perm);

  /* Apply vertex permutation to graph, optionally copying into user's result
   * array. */

  /* Get the size of the permutation fragment on each node. */
  count_type* perm_sizes = (count_type*)safe_malloc(size * sizeof(count_type));
  MPI_Allgather(&perm_local_size, 1, COUNT_MPI_TYPE, perm_sizes, 1, COUNT_MPI_TYPE, MPI_COMM_WORLD);

  /* Make a table consisting of the owner of each multiple of minpermsize (the
   * minimum permutation fragment size across all nodes).  Using this table,
   * only two locations need to be considered to find the owner of any given
   * vertex in the distributed permutation, avoiding a binary search.  Also,
   * sanity check the permutation to make sure its distribution is reasonable
   * (PRNG problems can cause bad distributions). */
  int minpermsize = perm_local_size;
  count_type i;
  for (i = 0; i < size; ++i) {
    if (perm_sizes[i] < minpermsize) minpermsize = perm_sizes[i];
  }

  if (minpermsize == 0) {
    fprintf(stderr, "One node has permutation size 0\n");
    abort();
  }
  count_type npermbuckets = (minpermsize == 0) ? N : (N + minpermsize - 1) / minpermsize;
  if (npermbuckets > size + 2 * N / perm_local_size) {
    fprintf(stderr, "Minimum permutation bucket size %d is too small compared to local size %" PRIcount_type "\n", minpermsize, perm_local_size);
    abort();
  }

  /* Prefix sum permutation sizes (exclusive) for owner lookups. */
  int* perm_displs = (int*)safe_malloc((size + 1) * sizeof(int));
  perm_displs[0] = 0;
  for (i = 1; i < size + 1; ++i) perm_displs[i] = perm_displs[i - 1] + perm_sizes[i - 1];
  count_type perm_total = perm_displs[size];

  /* Make the actual owner lookup table. */
  int* perm_owner_lookup = (int*)safe_malloc(npermbuckets * sizeof(int));
  int last_owner = 0;
  for (i = 0; i < npermbuckets; ++i) {
    count_type x = minpermsize * i;
    int j;
    for (j = last_owner; j < size; ++j) {
      if (x < perm_displs[j + 1] || j == size - 1) {
        perm_owner_lookup[i] = j;
        last_owner = j;
        break;
      }
    }
  }

  /* Macro to look up an element in the owner table. */
#define PERM_ELEMENT_OWNER(x) \
    (perm_owner_lookup[(x) / minpermsize] + ((x) >= perm_displs[perm_owner_lookup[(x) / minpermsize] + 1]))

  /* Create a bitmap (using characters) of which nodes need permutation
   * information from which other ones.  The used_dests map tells which nodes I
   * need permutation information from. */
  unsigned char* used_dests = (unsigned char*)safe_calloc(size, sizeof(unsigned char)); /* Uses zero init */
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  for (i = 0; i < *nedges; ++i) {
    if (local_edges[i].multiplicity == 0) continue;
    count_type x = local_edges[i].src;
    if (x >= perm_total) {fprintf(stderr, "x out of range: %" PRIcount_type "\n", x); abort();}
    used_dests[PERM_ELEMENT_OWNER(x)] = 1;
    x = local_edges[i].tgt;
    if (x >= perm_total) {fprintf(stderr, "x out of range: %" PRIcount_type "\n", x); abort();}
    used_dests[PERM_ELEMENT_OWNER(x)] = 1;
  }
#else
  for (i = 0; i < 2 * *nedges; ++i) {
    if (local_edges[i] == (count_type)(-1)) continue;
    count_type x = local_edges[i];
    if (x >= perm_total) {fprintf(stderr, "x out of range: %" PRIcount_type "\n", x); abort();}
    used_dests[PERM_ELEMENT_OWNER(x)] = 1;
  }
#endif
  /* The in_used_dests map tells which nodes I need to send permutation
   * information to. */
  unsigned char* in_used_dests = (unsigned char*)safe_malloc(size * sizeof(unsigned char));
  MPI_Alltoall(used_dests, 1, MPI_UNSIGNED_CHAR, in_used_dests, 1, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);

  /* Set up the MPI_Alltoallv inputs to multicast my permutation to all nodes
   * mentioned in in_used_dests and to receive from the nodes mentioned in
   * used_dests.  Also, set up the table containing the position of each
   * received permutation fragment in the temp array; this will be used to do
   * lookups later. */
  int* out_perm_counts = (int*)safe_malloc(size * sizeof(int));
  int* out_perm_displs = (int*)safe_calloc(size, sizeof(int)); /* Left at 0 */
  for (i = 0; i < size; ++i) {
    out_perm_counts[i] = in_used_dests[i] ? perm_local_size : 0;
  }
  int* in_perm_counts = (int*)safe_malloc(size * sizeof(int));
  int* in_perm_displs = (int*)safe_malloc((size + 1) * sizeof(int));
  in_perm_displs[0] = 0;
  for (i = 0; i < size; ++i) {
    in_perm_counts[i] = used_dests[i] ? perm_sizes[i] : 0;
    in_perm_displs[i + 1] = in_perm_displs[i] + in_perm_counts[i];
  }
  free(in_used_dests); in_used_dests = NULL;
  free(used_dests); used_dests = NULL;
  free(perm_sizes); perm_sizes = NULL;

  /* Move permutation data as necessary. */
  count_type* perm_temp = (count_type*)safe_malloc(in_perm_displs[size] * sizeof(count_type));
  MPI_Alltoallv(local_vertex_perm, out_perm_counts, out_perm_displs, COUNT_MPI_TYPE,
                perm_temp, in_perm_counts, in_perm_displs, COUNT_MPI_TYPE,
                MPI_COMM_WORLD);

  free(out_perm_counts); out_perm_counts = NULL;
  free(out_perm_displs); out_perm_displs = NULL;
  free(in_perm_counts); in_perm_counts = NULL;
  free(local_vertex_perm); local_vertex_perm = NULL;

  /* Permute endpoints of edges, copying into *result if necessary.  The
   * indexing into perm_temp here is to adjust from global offsets in the
   * permutation to offsets in an array containing the concatenation of the
   * fragments from only some nodes. */
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  *result = (count_type*)safe_malloc(2 * *nedges * sizeof(count_type));
  for (i = 0; i < *nedges; ++i) {
    if (local_edges[i].multiplicity != 0) {
      int owner = PERM_ELEMENT_OWNER(local_edges[i].src);
      (*result)[i * 2] = perm_temp[local_edges[i].src - perm_displs[owner] + in_perm_displs[owner]];
      owner = PERM_ELEMENT_OWNER(local_edges[i].tgt);
      (*result)[i * 2 + 1] = perm_temp[local_edges[i].tgt - perm_displs[owner] + in_perm_displs[owner]];
    } else {
      (*result)[i * 2] = (*result)[i * 2 + 1] = (count_type)(-1);
    }
  }
  free(local_edges); local_edges = NULL;
#else
  *result = local_edges;
  for (i = 0; i < 2 * *nedges; ++i) {
    if (local_edges[i] != (count_type)(-1)) {
      int owner = PERM_ELEMENT_OWNER(local_edges[i]);
      local_edges[i] = perm_temp[local_edges[i] - perm_displs[owner] + in_perm_displs[owner]];
    }
  }
#endif

  free(perm_temp); perm_temp = NULL;
  free(perm_displs); perm_displs = NULL;
  free(in_perm_displs); in_perm_displs = NULL;
  free(perm_owner_lookup); perm_owner_lookup = NULL;
#undef PERM_ELEMENT_OWNER
}
#endif

