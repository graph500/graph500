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
#include "oned_csr.h"
#include "onesided.h"
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

static oned_csr_graph g;
static unsigned long* g_in_queue;
static unsigned long* g_in_queue_summary;
static unsigned long* g_out_queue;
static unsigned long* g_out_queue_summary;
static unsigned long* g_visited;

const int ulong_bits = sizeof(unsigned long) * CHAR_BIT;
const int ulong_bits_squared = sizeof(unsigned long) * sizeof(unsigned long) * CHAR_BIT * CHAR_BIT;

static void allocate_memory(void) {
  int64_t maxlocalverts = g.max_nlocalverts;
  int64_t local_queue_summary_size = (maxlocalverts + ulong_bits_squared - 1) / ulong_bits_squared;
  int64_t local_queue_size = local_queue_summary_size * ulong_bits;
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
  convert_graph_to_oned_csr(tg, &g);
  allocate_memory(); /* Make sure all of the space is available */
  deallocate_memory();
}

void free_graph_data_structure(void) {
  free_oned_csr_graph(&g);
  /* deallocate_memory(); */
}

int bfs_writes_depth_map(void) {return 1;}

/* This version is the traditional level-synchronized BFS using two queues.  A
 * bitmap is used to indicate which vertices have been visited.  Messages are
 * sent and processed asynchronously throughout the code to hopefully overlap
 * communication with computation. */
void run_bfs(int64_t root, int64_t* pred) {
  allocate_memory();
  const ptrdiff_t nlocalverts = g.nlocalverts;
  const size_t* const restrict rowstarts = g.rowstarts;
  const int64_t* const restrict column = g.column;
  int64_t maxlocalverts = g.max_nlocalverts;

  /* Set up the visited bitmap. */
  const int ulong_bits = sizeof(unsigned long) * CHAR_BIT;
  const int ulong_bits_squared = ulong_bits * ulong_bits;
  int64_t local_queue_summary_size = (maxlocalverts + ulong_bits_squared - 1) / ulong_bits_squared;
  int64_t local_queue_size = local_queue_summary_size * ulong_bits;
  int lg_local_queue_size = lg_int64_t(local_queue_size);
  int64_t global_queue_summary_size = MUL_SIZE(local_queue_summary_size);
  int64_t global_queue_size = MUL_SIZE(local_queue_size);

#define SWIZZLE_VERTEX(c) ((VERTEX_OWNER(c) << lg_local_queue_size) * ulong_bits | VERTEX_LOCAL(c))
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

#define SET_IN(v) do {int64_t vs = SWIZZLE_VERTEX(v); size_t word_idx = vs / ulong_bits; int bit_idx = vs % ulong_bits; unsigned long mask = (1UL << bit_idx); in_queue_summary[word_idx / ulong_bits] |= (1UL << (word_idx % ulong_bits)); in_queue[word_idx] |= mask;} while (0)
#define TEST_IN(vs) (((in_queue_summary[vs / ulong_bits / ulong_bits] & (1UL << ((vs / ulong_bits) % ulong_bits))) != 0) && ((in_queue[vs / ulong_bits] & (1UL << (vs % ulong_bits))) != 0))
#define TEST_VISITED_LOCAL(v) ((visited[(v) / ulong_bits] & (1UL << ((v) % ulong_bits))) != 0)
// #define SET_VISITED_LOCAL(v) do {size_t word_idx = (v) / ulong_bits; int bit_idx = (v) % ulong_bits; unsigned long mask = (1UL << bit_idx); __sync_fetch_and_or(&visited[word_idx], mask); __sync_fetch_and_or(&out_queue[word_idx], mask);} while (0)
#define SET_VISITED_LOCAL(v) do {size_t word_idx = (v) / ulong_bits; int bit_idx = (v) % ulong_bits; unsigned long mask = (1UL << bit_idx); visited[word_idx] |= mask; out_queue[word_idx] |= mask;} while (0)

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
    ptrdiff_t i, ii;
#if 0
#pragma omp parallel for schedule(static)
    for (i = 0; i < global_queue_summary_size; ++i) {
      unsigned long val = 0UL;
      int j;
      unsigned long mask = 1UL;
      for (j = 0; j < ulong_bits; ++j, mask <<= 1) {
        if (in_queue[i * ulong_bits + j]) val |= mask;
      }
      in_queue_summary[i] = val;
    }
#endif
    unsigned long not_done = 0;
#pragma omp parallel for schedule(static) reduction(|:not_done)
    for (ii = 0; ii < nlocalverts; ii += ulong_bits) {
      size_t i, i_end = ii + ulong_bits;
      if (i_end > nlocalverts) i_end = nlocalverts;
      for (i = ii; i < i_end; ++i) {
        if (!TEST_VISITED_LOCAL(i)) {
          size_t j, j_end = rowstarts[i + 1];
          for (j = rowstarts[i]; j < j_end; ++j) {
            int64_t v1 = column[j];
            int64_t v1_swizzled = SWIZZLE_VERTEX(v1);
            if (TEST_IN(v1_swizzled)) {
              pred[i] = (v1 & INT64_C(0xFFFFFFFFFFFF)) | ((int64_t)cur_level << 48);
              not_done |= 1;
              SET_VISITED_LOCAL(i);
              break;
            }
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
      for (j = 0; j < ulong_bits; ++j, mask <<= 1) {
        unsigned long full_val = out_queue[i * ulong_bits + j];
        visited[i * ulong_bits + j] |= full_val;
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
    int64_t v = vertex[i];
    owner[i] = VERTEX_OWNER(v);
    local[i] = VERTEX_LOCAL(v);
  }
}

int64_t vertex_to_global_for_pred(int v_rank, size_t v_local) {
  return VERTEX_TO_GLOBAL(v_rank, v_local);
}

size_t get_nlocalverts_for_pred(void) {
  return g.nlocalverts;
}
