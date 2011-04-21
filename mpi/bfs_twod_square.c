/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#define _GNU_SOURCE
#include "common.h"
#include "twod_square.h"
#include "onesided.h"
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

static twod_square_graph g;
static unsigned long* g_in_queue;
static unsigned long* g_out_queue;
static unsigned long* g_visited;
static MPI_Comm comm_this_row, comm_this_column, comm_all_active;
static unsigned long* g_edge_flags; /* (Overapproximation of) whether each edge should become a pred map entry */

const int ulong_bits = sizeof(unsigned long) * CHAR_BIT;
const int ulong_bits_squared = sizeof(unsigned long) * sizeof(unsigned long) * CHAR_BIT * CHAR_BIT;

static void allocate_memory(void) {
  if (g.this_rank_active) {
    const int my_row = g.my_row;
    const int my_col = g.my_col;
    const int square_side = g.square_side;
    const size_t nlocalverts_1d = g.my_nlocalverts_1d;
    assert (nlocalverts_1d == get_nlocalverts_for_pred());
    const size_t nlocalverts_2d_src = g.my_nlocalverts_2d_src;
    const size_t nlocalverts_2d_tgt = g.my_nlocalverts_2d_tgt;
    const size_t nlocaledges = g.nlocaledges;

    int64_t local_in_queue_summary_size = (nlocalverts_2d_src + ulong_bits_squared - 1) / ulong_bits_squared;
    int64_t local_in_queue_size = local_in_queue_summary_size * ulong_bits;
    int64_t local_out_queue_summary_size = (nlocalverts_2d_tgt + ulong_bits_squared - 1) / ulong_bits_squared;
    int64_t local_out_queue_size = local_out_queue_summary_size * ulong_bits;
    g_in_queue = (unsigned long*)xmalloc(local_in_queue_size * sizeof(unsigned long)); /* By row */
    g_out_queue = (unsigned long*)xmalloc(local_out_queue_size * sizeof(unsigned long)); /* By column */
    g_visited = (unsigned long*)xmalloc(local_out_queue_size * sizeof(unsigned long)); /* By column to line up with out_queue */
    g_edge_flags = (unsigned long*)xmalloc((nlocaledges + ulong_bits - 1) / ulong_bits * sizeof(unsigned long));
  } else {
    g_in_queue = NULL;
    g_out_queue = NULL;
    g_visited = NULL;
    g_edge_flags = NULL;
  }
  MPI_Comm_split(MPI_COMM_WORLD, (g.this_rank_active ? g.my_row : MPI_UNDEFINED), g.my_col, &comm_this_row);
  MPI_Comm_split(MPI_COMM_WORLD, (g.this_rank_active ? g.my_col : MPI_UNDEFINED), g.my_row, &comm_this_column);
  MPI_Comm_split(MPI_COMM_WORLD, (g.this_rank_active ? 0 : MPI_UNDEFINED), 0, &comm_all_active);
}

static void deallocate_memory_1(void) { /* Skips communicators and edge flags */
  if (g_in_queue != NULL) {free(g_in_queue); g_in_queue = NULL;}
  if (g_out_queue != NULL) {free(g_out_queue); g_out_queue = NULL;}
  if (g_visited != NULL) {free(g_visited); g_visited = NULL;}
}

static void deallocate_memory_2(void) { /* Free rest of stuff */
  if (comm_this_row != MPI_COMM_NULL) MPI_Comm_free(&comm_this_row);
  if (comm_this_column != MPI_COMM_NULL) MPI_Comm_free(&comm_this_column);
  if (comm_all_active != MPI_COMM_NULL) MPI_Comm_free(&comm_all_active);
  if (g_edge_flags != NULL) {free(g_edge_flags); g_edge_flags = NULL;}
}

void make_graph_data_structure(const tuple_graph* const tg) {
  twod_square_graph g_temp;
  g = g_temp; /* Make the global variables uninitialized for valgrind */
  convert_graph_to_twod_square(tg, &g);
  allocate_memory(); /* Make sure all of the space is available */
  deallocate_memory_1();
  deallocate_memory_2();
}

void free_graph_data_structure(void) {
  free_twod_square_graph(&g);
  /* deallocate_memory(); */
}

void run_bfs(int64_t root, int64_t* pred) {
  allocate_memory();
  if (g.this_rank_active) {
    const int my_row = g.my_row;
    const int my_col = g.my_col;
    const int square_side = g.square_side;
    const size_t nlocalverts_1d = g.my_nlocalverts_1d;
    assert (nlocalverts_1d == get_nlocalverts_for_pred());
    const size_t nlocalverts_2d_src = g.my_nlocalverts_2d_src;
    const size_t nlocalverts_2d_tgt = g.my_nlocalverts_2d_tgt;
    const uint64_t* const restrict localedges = g.localedges;
    const size_t nlocaledges = g.nlocaledges;

    /* Set up the visited bitmap. */
    const int ulong_bits = ULONG_BITS;
    const int ulong_bits_squared = ULONG_BITS_SQUARED;

    int64_t local_in_queue_summary_size = nlocalverts_2d_src / ulong_bits_squared;
    int64_t local_in_queue_size = nlocalverts_2d_src / ulong_bits;
    int64_t local_out_queue_summary_size = nlocalverts_2d_tgt / ulong_bits_squared;
    int64_t local_out_queue_size = nlocalverts_2d_tgt / ulong_bits;
    unsigned long* restrict in_queue = g_in_queue;
    memset(in_queue, 0, local_in_queue_size * sizeof(unsigned long));
    unsigned long* restrict out_queue = g_out_queue;
    unsigned long* restrict visited = g_visited;
    memset(visited, 0, local_out_queue_size * sizeof(unsigned long));
    unsigned long* restrict edge_flags = g_edge_flags;
    memset(edge_flags, 0, (nlocaledges + ulong_bits - 1) / ulong_bits * sizeof(unsigned long));

#define SET_IN_LOCAL_SRC(v) do {size_t word_idx = (v) / ulong_bits; int bit_idx = (v) % ulong_bits; unsigned long mask = (1UL << bit_idx); in_queue[word_idx] |= mask;} while (0)
#define TEST_IN_LOCAL_SRC(v) ((in_queue[(v) / ulong_bits] & (1UL << ((v) % ulong_bits))) != 0)
#define TEST_VISITED_LOCAL_TGT(v) ((visited[(v) / ulong_bits] & (1UL << ((v) % ulong_bits))) != 0)
#define SET_VISITED_LOCAL_TGT(v) do {size_t word_idx = (v) / ulong_bits; int bit_idx = (v) % ulong_bits; unsigned long mask = (1UL << bit_idx); visited[word_idx] |= mask; out_queue[word_idx] |= mask;} while (0)
#define TAS_VISITED_LOCAL_TGT(v) (((__sync_fetch_and_or(&visited[(v) / ulong_bits], (1UL << ((v) % ulong_bits))) & (1UL << ((v) % ulong_bits))) != 0) ? 1 : (__sync_fetch_and_or(&out_queue[(v) / ulong_bits], (1UL << ((v) % ulong_bits))), 0))
#define SET_EDGE_FLAG(e) do {_Pragma("omp atomic") edge_flags[(e) / ulong_bits] |= (1UL << ((e) % ulong_bits));} while (0)

    if (vertex_owner_2d(root, &g) == g.my_row) {
      SET_IN_LOCAL_SRC(vertex_local_2d(root, &g));
    }
    if (vertex_owner_2d(root, &g) == g.my_col) {
      SET_VISITED_LOCAL_TGT(vertex_local_2d(root, &g));
    }
    {ptrdiff_t i; _Pragma("omp parallel for schedule(static)") for (i = 0; i < nlocalverts_1d; ++i) pred[i] = -1;}
    if (vertex_owner_1d(root, &g) == rank) {
      assert (vertex_to_global_1d(vertex_owner_1d(root, &g), vertex_local_1d(root, &g), &g) == root);
      assert (vertex_local_1d(root, &g) < nlocalverts_1d);
      pred[vertex_local_1d(root, &g)] = root;
    }
    uint16_t cur_level = 0;
    while (1) {
      ++cur_level;
      if (rank == 0) fprintf(stderr, "BFS level %" PRIu16 "\n", cur_level);
      memset(out_queue, 0, local_out_queue_size * sizeof(unsigned long));
      ptrdiff_t i, ii;
      unsigned long not_done = 0;
#pragma omp parallel for schedule(static) reduction(|:not_done)
      for (ii = 0; ii < nlocaledges; ++ii) {
        uint64_t e = localedges[ii];
        uint32_t local_src = (uint32_t)(e >> 32);
        uint32_t local_tgt = (uint32_t)(e & UINT32_C(0xFFFFFFFF));
        assert (local_src < nlocalverts_2d_src);
        assert (local_tgt < nlocalverts_2d_tgt);
        if (TEST_IN_LOCAL_SRC(local_src) && !TAS_VISITED_LOCAL_TGT(local_tgt)) {
          SET_EDGE_FLAG(ii);
          not_done |= 1;
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &not_done, 1, MPI_UNSIGNED_LONG, MPI_BOR, comm_all_active);
      if (not_done == 0) break;
      MPI_Allreduce(MPI_IN_PLACE, out_queue, local_out_queue_size, MPI_UNSIGNED_LONG, MPI_BOR, comm_this_column);
#pragma omp parallel for
      for (i = 0; i < local_out_queue_size; ++i) {
        visited[i] |= out_queue[i];
      }
      if (my_row == my_col) {
        assert (nlocalverts_2d_src == nlocalverts_2d_tgt);
        assert (local_in_queue_size == local_out_queue_size);
        memcpy(in_queue, out_queue, local_in_queue_size * sizeof(unsigned long));
      }
      MPI_Bcast(in_queue, local_in_queue_size, MPI_UNSIGNED_LONG, my_row, comm_this_row);
    }
    deallocate_memory_1();
    int* edge_flags_per_dest = (int*)xmalloc(square_side * square_side * sizeof(int));
    int* edge_flag_displs_per_dest = (int*)xmalloc((square_side * square_side + 1) * sizeof(int));
    memset(edge_flags_per_dest, 0, square_side * square_side * sizeof(int));
    edge_flag_displs_per_dest[0] = 0;
    size_t ii;
#pragma omp parallel for
    for (ii = 0; ii < nlocaledges; ++ii) {
      if (edge_flags[ii / ulong_bits] & (1UL << (ii % ulong_bits))) {
        uint64_t e = localedges[ii];
        uint32_t local_tgt = (uint32_t)(e & UINT32_C(0xFFFFFFFF));
#pragma omp atomic
        ++edge_flags_per_dest[vertex_owner_1d(vertex_to_global_2d(my_col, local_tgt, &g), &g)];
      }
    }
    for (ii = 0; ii < square_side * square_side; ++ii) {
      edge_flag_displs_per_dest[ii + 1] = edge_flag_displs_per_dest[ii] + edge_flags_per_dest[ii];
    }
    uint64_t* outdata = (uint64_t*)xmalloc((size_t)edge_flag_displs_per_dest[square_side * square_side] * sizeof(uint64_t));
    int* edge_flags_per_src = (int*)xmalloc(square_side * square_side * sizeof(int));
    MPI_Alltoall(edge_flags_per_dest, 1, MPI_INT, edge_flags_per_src, 1, MPI_INT, comm_all_active);
    int* edge_flag_displs_per_src = (int*)xmalloc((square_side * square_side + 1) * sizeof(int));
    edge_flag_displs_per_src[0] = 0;
    for (ii = 0; ii < square_side * square_side; ++ii) {
      edge_flag_displs_per_src[ii + 1] = edge_flag_displs_per_src[ii] + edge_flags_per_src[ii];
    }
    int* edge_flag_inserts_per_dest = (int*)xmalloc(square_side * square_side * sizeof(int));
    memcpy(edge_flag_inserts_per_dest, edge_flag_displs_per_dest, square_side * square_side * sizeof(int));
#pragma omp parallel for
    for (ii = 0; ii < nlocaledges; ++ii) {
      if (edge_flags[ii / ulong_bits] & (1UL << (ii % ulong_bits))) {
        uint64_t e = localedges[ii];
        uint32_t local_tgt = (uint32_t)(e & UINT32_C(0xFFFFFFFF));
        outdata[__sync_fetch_and_add(&edge_flag_inserts_per_dest[vertex_owner_1d(vertex_to_global_2d(my_col, local_tgt, &g), &g)], 1)] = e;
      }
    }
    free(edge_flag_inserts_per_dest); edge_flag_inserts_per_dest = NULL;
    uint64_t* indata = (uint64_t*)xmalloc((size_t)edge_flag_displs_per_src[square_side * square_side] * sizeof(uint64_t));
    MPI_Alltoallv(outdata, edge_flags_per_dest, edge_flag_displs_per_dest, MPI_UINT64_T,
                  indata, edge_flags_per_src, edge_flag_displs_per_src, MPI_UINT64_T,
                  comm_all_active);
    free(outdata); outdata = NULL;
    free(edge_flags_per_dest); edge_flags_per_dest = NULL;
    free(edge_flag_displs_per_dest); edge_flag_displs_per_dest = NULL;
    free(edge_flags_per_src); edge_flags_per_src = NULL;
    free(edge_flag_displs_per_src); edge_flag_displs_per_src = NULL;
#pragma omp parallel for schedule(dynamic)
    for (ii = 0; ii < square_side * square_side; ++ii) {
      int i = edge_flag_displs_per_src[ii], i_end = edge_flag_displs_per_src[ii + 1];
      for (; i < i_end; ++i) {
        uint64_t e = indata[i];
        uint32_t their_local_src = (uint32_t)(e >> 32);
        uint32_t their_local_tgt = (uint32_t)(e & UINT32_C(0xFFFFFFFF));
        int their_row = ii / square_side;
        int their_col = ii % square_side;
        uint32_t my_local_tgt = vertex_local_1d(vertex_to_global_2d(their_col, their_local_tgt, &g), &g);
        pred[my_local_tgt] = vertex_to_global_2d(their_row, their_local_src, &g);
      }
    }
    free(indata); indata = NULL;
  }
  deallocate_memory_2();
}

void get_vertex_distribution_for_pred(size_t count, const int64_t* vertex_p, int* owner_p, size_t* local_p) {
  const int64_t* restrict vertex = vertex_p;
  int* restrict owner = owner_p;
  size_t* restrict local = local_p;
  ptrdiff_t i;
#pragma omp parallel for
  for (i = 0; i < (ptrdiff_t)count; ++i) {
    int64_t v = vertex[i];
    if (v < 0 || v > (INT64_C(1) << g.lg_nglobalverts)) {
      owner[i] = MPI_PROC_NULL;
      local[i] = -1;
    } else {
      owner[i] = vertex_owner_1d(v, &g);
      local[i] = vertex_local_1d(v, &g);
      assert (owner[i] >= 0 && owner[i] < g.square_side * g.square_side);
      assert (local[i] < vertex_offset_for_coord_1d(owner[i] + 1, &g) - vertex_offset_for_coord_1d(owner[i], &g));
    }
  }
}

int64_t vertex_to_global_for_pred(int v_rank, size_t v_local) {
  int64_t v = vertex_to_global_1d(v_rank, v_local, &g);
  assert (v >= 0 && v < (INT64_C(1) << g.lg_nglobalverts));
  assert (vertex_owner_1d(v, &g) == v_rank);
  assert (vertex_local_1d(v, &g) == v_local);
  return v;
}

size_t get_nlocalverts_for_pred(void) {
  if (g.this_rank_active) assert (vertex_offset_for_coord_1d(rank + 1, &g) - vertex_offset_for_coord_1d(rank, &g) == g.my_nlocalverts_1d);
  return g.this_rank_active ? g.my_nlocalverts_1d : 0;
}
