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
#include "apply_permutation_mpi.h"
#include "scramble_edges.h"
#include "utils.h"
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

/* This version is for sequential machines, OpenMP, and the XMT. */
void scramble_edges_shared(uint64_t userseed1, uint64_t userseed2, int64_t nedges, int64_t* result /* Input and output array of edges (size = 2 * nedges) */) {
  mrg_state st;
  uint_fast32_t seed[5];
  int64_t* perm = (int64_t*)xmalloc(nedges * sizeof(int64_t));
  make_mrg_seed(userseed1, userseed2, seed);
  mrg_seed(&st, seed);
  mrg_skip(&st, 5, 0, 0); /* To make offset different from other PRNG uses */
  rand_sort_shared(&st, nedges, perm);
  int64_t* new_result = (int64_t*)xmalloc(nedges * 2 * sizeof(int64_t));
  int64_t i;
#ifdef __MTA__
#pragma mta assert parallel
#pragma mta block schedule
#endif
#ifdef GRAPH_GENERATOR_OMP
#pragma omp parallel for
#endif
  for (i = 0; i < nedges; ++i) {
    int64_t p = perm[i];
    new_result[i * 2 + 0] = result[p * 2 + 0];
    new_result[i * 2 + 1] = result[p * 2 + 1];
  }
  free(perm);
  memcpy(result, new_result, nedges * 2 * sizeof(int64_t));
  free(new_result);
}

#ifdef GRAPH_GENERATOR_MPI
/* For MPI distributed memory. */
void scramble_edges_mpi(MPI_Comm comm,
                        const uint64_t userseed1, const uint64_t userseed2,
                        const int64_t local_nedges_in,
                        const int64_t* const local_edges_in,
                        int64_t* const local_nedges_out_ptr,
                        int64_t** const local_edges_out_ptr /* Allocated using xmalloc() by scramble_edges_mpi */) {
  int rank, size;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  mrg_state st;
  uint_fast32_t seed[5];
  make_mrg_seed(userseed1, userseed2, seed);
  mrg_seed(&st, seed);
  mrg_skip(&st, 5, 0, 0); /* To make offset different from other PRNG uses */
  int64_t total_nedges;
  MPI_Allreduce((void*)&local_nedges_in, &total_nedges, 1, INT64_T_MPI_TYPE, MPI_SUM, comm);
  int64_t local_nedges_out; /* = local permutation size */
  int64_t* local_perm;
  rand_sort_mpi(comm, &st, total_nedges, &local_nedges_out, &local_perm);
  *local_nedges_out_ptr = local_nedges_out;

  /* Gather permutation information and fast owner lookup cache (code in
   * apply_permutation_mpi.c). */
  int64_t* edge_displs = (int64_t*)xmalloc((size + 1) * sizeof(int64_t));
  int* edge_owner_table;
  int64_t* edge_owner_cutoff;
  int lg_minedgecount;
  int64_t maxedgecount;
  gather_block_distribution_info(comm, local_nedges_in, total_nedges, edge_displs, &edge_owner_table, &edge_owner_cutoff, &lg_minedgecount, &maxedgecount);

  /* Originally from apply_permutation_mpi.c */
#define LOOKUP_EDGE_OWNER(v) \
  (edge_owner_table[(v) >> lg_minedgecount] + \
   ((v) >= edge_owner_cutoff[(v) >> lg_minedgecount]))

  /* Apply permutation.  Output distribution is same as distribution of
   * generated edge permutation. */

  /* Count number of requests to send to each destination. */
  int* send_counts = (int*)xcalloc(size, sizeof(int)); /* Uses zero-init */
  int64_t i;
  for (i = 0; i < local_nedges_out; ++i) {
    ++send_counts[LOOKUP_EDGE_OWNER(local_perm[i])];
  }

  /* Prefix sum to get displacements. */
  int* send_displs = (int*)xmalloc((size + 1) * sizeof(int));
  send_displs[0] = 0;
  for (i = 0; i < size; ++i) {
    send_displs[i + 1] = send_displs[i] + send_counts[i];
  }
  assert (send_displs[size] == local_nedges_out);

  /* Put edges into buffer by destination; also keep around index values for
   * where to write the result. */
  int64_t* sendbuf = (int64_t*)xmalloc(local_nedges_out * sizeof(int64_t));
  int64_t* reply_loc_buf = (int64_t*)xmalloc(local_nedges_out * sizeof(int64_t));
  int* send_offsets = (int*)xmalloc((size + 1) * sizeof(int));
  memcpy(send_offsets, send_displs, (size + 1) * sizeof(int));
  for (i = 0; i < local_nedges_out; ++i) {
    int write_index = send_offsets[LOOKUP_EDGE_OWNER(local_perm[i])];
    sendbuf[write_index] = local_perm[i];
    reply_loc_buf[write_index] = i;
    ++send_offsets[LOOKUP_EDGE_OWNER(local_perm[i])];
  }
  for (i = 0; i < size; ++i) assert (send_offsets[i] == send_displs[i + 1]);
  free(send_offsets); send_offsets = NULL;
  free(local_perm); local_perm = NULL;
#undef LOOKUP_EDGE_OWNER
  free(edge_owner_table); edge_owner_table = NULL;
  free(edge_owner_cutoff); edge_owner_cutoff = NULL;

  /* Find out how many requests I will be receiving. */
  int* recv_counts = (int*)xmalloc(size * sizeof(int));
  MPI_Alltoall(send_counts, 1, MPI_INT,
               recv_counts, 1, MPI_INT,
               comm);

  /* Compute their displacements. */
  int* recv_displs = (int*)xmalloc((size + 1) * sizeof(int));
  recv_displs[0] = 0;
  for (i = 0; i < size; ++i) {
    recv_displs[i + 1] = recv_displs[i] + recv_counts[i];
  }

  /* Make receive and reply buffers. */
  int64_t* recvbuf = (int64_t*)xmalloc(recv_displs[size] * sizeof(int64_t));
  int64_t* replybuf = (int64_t*)xmalloc(recv_displs[size] * 2 * sizeof(int64_t));

  /* Move requests for edges into receive buffer. */
  MPI_Alltoallv(sendbuf, send_counts, send_displs, INT64_T_MPI_TYPE,
                recvbuf, recv_counts, recv_displs, INT64_T_MPI_TYPE,
                comm);
  free(sendbuf); sendbuf = NULL;

  /* Put requested edges into response buffer. */
  int64_t my_edge_offset = edge_displs[rank];
  for (i = 0; i < recv_displs[size]; ++i) {
    replybuf[i * 2 + 0] = local_edges_in[(recvbuf[i] - my_edge_offset) * 2 + 0];
    replybuf[i * 2 + 1] = local_edges_in[(recvbuf[i] - my_edge_offset) * 2 + 1];
  }
  free(recvbuf); recvbuf = NULL;
  free(edge_displs); edge_displs = NULL;

  /* Send replies back. */
  int64_t* reply_edges = (int64_t*)xmalloc(local_nedges_out * 2 * sizeof(int64_t));
  for (i = 0; i < size; ++i) { /* Sending back two values for each request */
    recv_counts[i] *= 2;
    recv_displs[i] *= 2;
    send_counts[i] *= 2;
    send_displs[i] *= 2;
  }
  MPI_Alltoallv(replybuf, recv_counts, recv_displs, INT64_T_MPI_TYPE,
                reply_edges, send_counts, send_displs, INT64_T_MPI_TYPE,
                comm);
  free(replybuf); replybuf = NULL;
  free(recv_counts); recv_counts = NULL;
  free(recv_displs); recv_displs = NULL;
  free(send_counts); send_counts = NULL;
  free(send_displs); send_displs = NULL;

  /* Make output array of edges. */
  int64_t* local_edges_out = (int64_t*)xmalloc(local_nedges_out * 2 * sizeof(int64_t));
  *local_edges_out_ptr = local_edges_out;

  /* Put edges into output array. */
  for (i = 0; i < local_nedges_out; ++i) {
    local_edges_out[reply_loc_buf[i] * 2 + 0] = reply_edges[2 * i + 0];
    local_edges_out[reply_loc_buf[i] * 2 + 1] = reply_edges[2 * i + 1];
  }
  free(reply_loc_buf); reply_loc_buf = NULL;
  free(reply_edges); reply_edges = NULL;
}
#endif
