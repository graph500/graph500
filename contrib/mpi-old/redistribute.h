/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef REDISTRIBUTE_H
#define REDISTRIBUTE_H

/* Macro that is the equivalent of a C++ template to implement generic edge
 * redistribution. */

#include "common.h"
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#define MAKE_REDISTRIBUTE_FUNC_COUNT_EDGE_BY_OWNER(RF_OWNER, RF_WRITE) \
  _Pragma("omp atomic") ++edge_counts_per_owner[(RF_OWNER)];

#define MAKE_REDISTRIBUTE_FUNC_STORE_EDGE_TO_BUFFER(RF_OWNER, RF_WRITE) \
  { \
    ptrdiff_t pos = __sync_fetch_and_add(&edge_inserts_per_owner[(RF_OWNER)], 1); \
    RF_WRITE((&edges_to_send[pos]), (v0), (v1)); \
  }

#define MAKE_REDISTRIBUTE_FUNC(RF_FUNCNAME, RF_EXTRA_PARAMS, RF_DECLARE_AND_INIT_GRAPH_SO_FAR, RF_CALL_ON_EDGES, RF_EDGE_BUFFER_TYPE, RF_EDGE_BUFFER_MPI_TYPE, RF_PRECOMPRESS_INCOMING_DATA, RF_MERGE_INTO_GRAPH_SO_FAR, RF_FREE_PRECOMPRESSED_DATA, RF_BUILD_FINAL_DATA_STRUCTURE_FROM_GRAPH_SO_FAR, RF_CLEAR_GRAPH_SO_FAR) \
void RF_FUNCNAME(const tuple_graph* const tg, RF_EXTRA_PARAMS) { \
  /* Scan for vertex count and edges per destination rank. */ \
  uint64_t max_vertex = UINT64_C(0); \
  int* edge_counts_per_owner = (int*)xmalloc(size * sizeof(int)); \
  int* edge_counts_per_sender = (int*)xmalloc(size * sizeof(int)); \
  int* edge_displs_per_owner = (int*)xmalloc((size + 1) * sizeof(int)); \
  int* edge_displs_per_sender = (int*)xmalloc((size + 1) * sizeof(int)); \
  int* edge_inserts_per_owner = (int*)xmalloc(size * sizeof(int)); \
  size_t block_count = ITERATE_TUPLE_GRAPH_BLOCK_COUNT(tg); \
  RF_DECLARE_AND_INIT_GRAPH_SO_FAR \
  int restart = 0; \
  do { \
    restart = 0; \
    ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize) { \
      ptrdiff_t i; \
      uint64_t prev_max_vertex = max_vertex; \
      /* The |= trick (to allow an OMP reduction) works because we round max_vertex \
       * up to one less than a power of 2. */ \
      _Pragma("omp parallel for reduction(|:max_vertex)") \
      for (i = 0; i < bufsize; ++i) { \
        int64_t v0 = get_v0_from_edge(&buf[i]); \
        int64_t v1 = get_v1_from_edge(&buf[i]); \
        if (v0 == v1) continue; \
        max_vertex |= (uint64_t)(v0 | v1); \
      } \
      MPI_Allreduce(MPI_IN_PLACE, &max_vertex, 1, MPI_UINT64_T, MPI_BOR, MPI_COMM_WORLD); \
      while (/* max_vertex not of the form 2**k-1 */ (max_vertex & (max_vertex + 1)) != 0) { \
        /* Set the lowest unset bit in max_vertex. */ \
        max_vertex |= max_vertex + 1; \
      } \
      if (prev_max_vertex != 0 && max_vertex != prev_max_vertex) { \
        if (rank == 0) { \
          fprintf(stderr, "Restarting because of change of max_vertex from %" PRIu64 " to %" PRIu64 "\n", prev_max_vertex, max_vertex); \
        } \
        restart = 1; \
        RF_CLEAR_GRAPH_SO_FAR \
        ITERATE_TUPLE_GRAPH_BREAK; \
      } \
      int lg_nglobalverts = 0; \
      for (; (max_vertex >> lg_nglobalverts) != 0; ++lg_nglobalverts) {} \
      memset(edge_counts_per_owner, 0, size * sizeof(int)); \
      _Pragma("omp parallel for") \
      for (i = 0; i < bufsize; ++i) { \
        int64_t v0 = get_v0_from_edge(&buf[i]); \
        int64_t v1 = get_v1_from_edge(&buf[i]); \
        if (v0 == v1) continue; \
        RF_CALL_ON_EDGES(v0, v1, lg_nglobalverts, MAKE_REDISTRIBUTE_FUNC_COUNT_EDGE_BY_OWNER) \
      } \
      MPI_Alltoall(edge_counts_per_owner, 1, MPI_INT, edge_counts_per_sender, 1, MPI_INT, MPI_COMM_WORLD); \
      edge_displs_per_owner[0] = 0; \
      edge_displs_per_sender[0] = 0; \
      for (i = 0; i < size; ++i) { \
        edge_displs_per_owner[i + 1] = edge_displs_per_owner[i] + edge_counts_per_owner[i]; \
        edge_displs_per_sender[i + 1] = edge_displs_per_sender[i] + edge_counts_per_sender[i]; \
      } \
      memcpy(edge_inserts_per_owner, edge_displs_per_owner, size * sizeof(int)); \
      RF_EDGE_BUFFER_TYPE* edges_to_send = (RF_EDGE_BUFFER_TYPE*)xMPI_Alloc_mem((size_t)edge_displs_per_owner[size] * sizeof(RF_EDGE_BUFFER_TYPE)); \
      _Pragma("omp parallel for") \
      for (i = 0; i < bufsize; ++i) { \
        int64_t v0 = get_v0_from_edge(&buf[i]); \
        int64_t v1 = get_v1_from_edge(&buf[i]); \
        if (v0 == v1) continue; \
        RF_CALL_ON_EDGES(v0, v1, lg_nglobalverts, MAKE_REDISTRIBUTE_FUNC_STORE_EDGE_TO_BUFFER) \
      } \
      ITERATE_TUPLE_GRAPH_RELEASE_BUFFER; \
      size_t edges_received_this_block = (size_t)edge_displs_per_sender[size]; \
      RF_EDGE_BUFFER_TYPE* restrict edges_to_recv = (RF_EDGE_BUFFER_TYPE*)xMPI_Alloc_mem(edges_received_this_block * sizeof(RF_EDGE_BUFFER_TYPE)); \
   \
      MPI_Alltoallv(edges_to_send, edge_counts_per_owner, edge_displs_per_owner, RF_EDGE_BUFFER_MPI_TYPE, \
                    edges_to_recv, edge_counts_per_sender, edge_displs_per_sender, RF_EDGE_BUFFER_MPI_TYPE, \
                    MPI_COMM_WORLD); \
      MPI_Free_mem(edges_to_send); \
   \
      { \
        RF_PRECOMPRESS_INCOMING_DATA(lg_nglobalverts, edges_to_recv, edges_received_this_block) \
        MPI_Free_mem(edges_to_recv); \
        RF_MERGE_INTO_GRAPH_SO_FAR \
        RF_FREE_PRECOMPRESSED_DATA \
      } \
    } ITERATE_TUPLE_GRAPH_END; \
  } while (restart); \
  free(edge_counts_per_owner); \
  free(edge_counts_per_sender); \
  free(edge_displs_per_owner); \
  free(edge_displs_per_sender); \
  free(edge_inserts_per_owner); \
 \
  RF_BUILD_FINAL_DATA_STRUCTURE_FROM_GRAPH_SO_FAR /* Should free graph_so_far as well */ \
}

#endif /* REDISTRIBUTE_H */
