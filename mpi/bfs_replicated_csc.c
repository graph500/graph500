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
#include "oned_csc.h"
#include "onesided.h"
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

static oned_csc_graph g;
static unsigned long* g_in_queue;
static unsigned long* g_in_queue_summary;
static unsigned long* g_out_queue;
static unsigned long* g_out_queue_summary;
static unsigned long* g_visited;

static void allocate_memory(void) {
  int64_t maxlocalverts = g.max_nlocalverts;
  int64_t local_queue_summary_size = (maxlocalverts + ULONG_BITS * ULONG_BITS - 1) / ULONG_BITS / ULONG_BITS;
  int64_t local_queue_size = local_queue_summary_size * ULONG_BITS;
  int64_t global_queue_summary_size = MUL_SIZE(local_queue_summary_size);
  int64_t global_queue_size = MUL_SIZE(local_queue_size);
  g_in_queue = (unsigned long*)xmalloc(global_queue_size * sizeof(unsigned long));
  g_in_queue_summary = (unsigned long*)xmalloc(global_queue_summary_size * sizeof(unsigned long));
  g_out_queue = (unsigned long*)xmalloc(local_queue_size * sizeof(unsigned long));
  g_out_queue_summary = (unsigned long*)xmalloc(local_queue_summary_size * sizeof(unsigned long));
  g_visited = (unsigned long*)xmalloc(local_queue_size * sizeof(unsigned long));
}

static void deallocate_memory(void) {
  free(g_in_queue); g_in_queue = NULL;
  free(g_in_queue_summary); g_in_queue_summary = NULL;
  free(g_out_queue); g_out_queue = NULL;
  free(g_out_queue_summary); g_out_queue_summary = NULL;
  free(g_visited); g_visited = NULL;
}

void make_graph_data_structure(const tuple_graph* const tg) {
  convert_graph_to_oned_csc(tg, &g);
  allocate_memory(); /* Make sure all of the space is available */
  deallocate_memory();
}

void free_graph_data_structure(void) {
  free_oned_csc_graph(&g);
  /* deallocate_memory(); */
}

int bfs_writes_depth_map(void) {
  return 0;
}

/* This version is the traditional level-synchronized BFS using two queues.  A
 * bitmap is used to indicate which vertices have been visited.  Messages are
 * sent and processed asynchronously throughout the code to hopefully overlap
 * communication with computation. */
void run_bfs(int64_t root, int64_t* pred) {
  allocate_memory();
  const ptrdiff_t nlocalverts = g.nlocalverts;
  /* const int64_t nglobalverts = g.nglobalverts; */
  const size_t* const restrict rowstarts = g.rowstarts;
  const int64_t* const restrict column = g.column;

  /* Set up the visited bitmap. */
  int lg_local_queue_size = g.lg_local_queue_size;
  int64_t local_queue_size = INT64_C(1) << lg_local_queue_size;
  int64_t local_queue_summary_size = local_queue_size / ULONG_BITS;
  int64_t global_queue_summary_size = MUL_SIZE(local_queue_summary_size);
  int64_t global_queue_size = MUL_SIZE(local_queue_size);

#if 0
  int64_t* restrict column_swizzled = (int64_t*)xmalloc(nlocaledges * sizeof(int64_t));
  {
    size_t i;
    for (i = 0; i < nlocaledges; ++i) {
      int64_t c = column[i];
      column_swizzled[i] = SWIZZLE_VERTEX(c);
    }
  }
#endif

  unsigned long* restrict in_queue = g_in_queue;
  memset(in_queue, 0, global_queue_size * sizeof(unsigned long));
  unsigned long* restrict in_queue_summary = g_in_queue_summary;
  memset(in_queue_summary, 0, global_queue_summary_size * sizeof(unsigned long));
  unsigned long* restrict out_queue = g_out_queue;
  unsigned long* restrict out_queue_summary = g_out_queue_summary;
  unsigned long* restrict visited = g_visited;
  memset(visited, 0, local_queue_size * sizeof(unsigned long));

#define SET_IN(v) do {int64_t vs = SWIZZLE_VERTEX(v); size_t word_idx = vs / ULONG_BITS; int bit_idx = vs % ULONG_BITS; unsigned long mask = (1UL << bit_idx); in_queue_summary[word_idx / ULONG_BITS] |= (1UL << (word_idx % ULONG_BITS)); in_queue[word_idx] |= mask;} while (0)
#define TEST_IN(vs) (((in_queue_summary[vs / ULONG_BITS / ULONG_BITS] & (1UL << ((vs / ULONG_BITS) % ULONG_BITS))) != 0) && ((in_queue[vs / ULONG_BITS] & (1UL << (vs % ULONG_BITS))) != 0))
#define TEST_VISITED_LOCAL(v) ((visited[(v) / ULONG_BITS] & (1UL << ((v) % ULONG_BITS))) != 0)
#define TAS_VISITED_LOCAL(v) (((__sync_fetch_and_or(&visited[(v) / ULONG_BITS], (1UL << ((v) % ULONG_BITS))) & (1UL << ((v) % ULONG_BITS))) != 0) ? 1 : (__sync_fetch_and_or(&out_queue[(v) / ULONG_BITS], (1UL << ((v) % ULONG_BITS))), 0))
// #define SET_VISITED_LOCAL(v) do {size_t word_idx = (v) / ULONG_BITS; int bit_idx = (v) % ULONG_BITS; unsigned long mask = (1UL << bit_idx); __sync_fetch_and_or(&visited[word_idx], mask); __sync_fetch_and_or(&out_queue[word_idx], mask);} while (0)
#define SET_VISITED_LOCAL(v) do {size_t word_idx = (v) / ULONG_BITS; int bit_idx = (v) % ULONG_BITS; unsigned long mask = (1UL << bit_idx); visited[word_idx] |= mask; out_queue[word_idx] |= mask;} while (0)

  SET_IN(root);
  {ptrdiff_t i; _Pragma("omp parallel for schedule(static)") for (i = 0; i < nlocalverts; ++i) pred[i] = -1;}
  if (VERTEX_OWNER(root) == rank) {
    pred[VERTEX_LOCAL(root)] = root;
    SET_VISITED_LOCAL(VERTEX_LOCAL(root));
  }
  uint16_t cur_level = 0;
  while (1) {
    ++cur_level;
#if 0
    if (rank == 0) fprintf(stderr, "BFS level %" PRIu16 "\n", cur_level);
#endif
    memset(out_queue, 0, local_queue_size * sizeof(unsigned long));
    // memset(out_queue_summary, 0, local_queue_summary_size * sizeof(unsigned long));
    ptrdiff_t i, ii_summary;
#if 0
#pragma omp parallel for schedule(static)
    for (i = 0; i < global_queue_summary_size; ++i) {
      unsigned long val = 0UL;
      int j;
      unsigned long mask = 1UL;
      for (j = 0; j < ULONG_BITS; ++j, mask <<= 1) {
        if (in_queue[i * ULONG_BITS + j]) val |= mask;
      }
      in_queue_summary[i] = val;
    }
#endif
    unsigned long not_done = 0;
#pragma omp parallel for schedule(static) reduction(|:not_done)
    for (ii_summary = 0; ii_summary < global_queue_summary_size; ++ii_summary) {
      uint64_t val_summary = in_queue_summary[ii_summary];
      if (val_summary == 0) continue;
      int ii_offset;
      ptrdiff_t ii;
      for (ii_offset = 0; ii_offset < ULONG_BITS; ++ii_offset) {
        if ((val_summary & (UINT64_C(1) << ii_offset)) == 0) continue;
        ii = ii_summary * ULONG_BITS + ii_offset;
        uint64_t val = in_queue[ii];
        if (val == 0) continue;
        size_t i, i_end = rowstarts[ii + 1];
        for (i = rowstarts[ii]; i < i_end; ++i) {
          int64_t c = column[i];
          int64_t v0_local = c / ULONG_BITS;
          if ((val & (UINT64_C(1) << (c % ULONG_BITS))) != 0 /* TEST_IN(v1_swizzled) */ && !TAS_VISITED_LOCAL(v0_local)) {
            assert (pred[v0_local] == -1);
            int64_t v1_swizzled = (int64_t)ii * ULONG_BITS + c % ULONG_BITS;
            pred[v0_local] = UNSWIZZLE_VERTEX(v1_swizzled);
            not_done |= 1;
          }
        }
      }
    }
#if 1
#pragma omp parallel for schedule(static)
    for (i = 0; i < local_queue_summary_size; ++i) {
      unsigned long val = 0UL;
      int j;
      unsigned long mask = 1UL;
      for (j = 0; j < ULONG_BITS; ++j, mask <<= 1) {
        unsigned long full_val = out_queue[i * ULONG_BITS + j];
        visited[i * ULONG_BITS + j] |= full_val;
        if (full_val) val |= mask;
      }
      out_queue_summary[i] = val;
      // not_done |= val;
    }
#endif
    MPI_Allreduce(MPI_IN_PLACE, &not_done, 1, MPI_UNSIGNED_LONG, MPI_BOR, MPI_COMM_WORLD);
    if (not_done == 0) break;
    MPI_Allgather(out_queue, local_queue_size, MPI_UNSIGNED_LONG, in_queue, local_queue_size, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MPI_Allgather(out_queue_summary, local_queue_summary_size, MPI_UNSIGNED_LONG, in_queue_summary, local_queue_summary_size, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  }
  deallocate_memory();
}

void get_vertex_distribution_for_pred(size_t count, const int64_t* vertex_p, int* owner_p, size_t* local_p) {
  const int64_t* restrict vertex = vertex_p;
  int* restrict owner = owner_p;
  size_t* restrict local = local_p;
  ptrdiff_t i;
#pragma omp parallel for
  for (i = 0; i < (ptrdiff_t)count; ++i) {
    owner[i] = VERTEX_OWNER(vertex[i]);
    local[i] = VERTEX_LOCAL(vertex[i]);
  }
}

int64_t vertex_to_global_for_pred(int v_rank, size_t v_local) {
  return VERTEX_TO_GLOBAL(v_rank, v_local);
}

size_t get_nlocalverts_for_pred(void) {
  return g.nlocalverts;
}
