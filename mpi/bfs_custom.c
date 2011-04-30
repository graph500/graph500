/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include "common.h"
#include "oned_csr.h"
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

/* Add your own BFS code into this file (or a copy of it). */

/* Data structure definitions: customize these for your own data distribution
 * and temporary data structures. */
static oned_csr_graph g;

void make_graph_data_structure(const tuple_graph* const tg) {
  convert_graph_to_oned_csr(tg, &g);
}

void free_graph_data_structure(void) {
  free_oned_csr_graph(&g);
}

int bfs_writes_depth_map(void) {
  /* Change to 1 if high 16 bits of each entry of pred are the (zero-based) BFS
   * level number, with UINT16_MAX for unreachable vertices. */
  return 0;
}

/* BFS implementation. */
void run_bfs(int64_t root, int64_t* pred) {
  /* Predefined entities you can use in your BFS (from common.h and oned_csr.h):
   *   + rank: global variable containing MPI rank
   *   + size: global variable containing MPI size
   *   + DIV_SIZE: single-parameter macro that divides by size (using a shift
   *     when properly set up)
   *   + MOD_SIZE: single-parameter macro that reduces modulo size (using a
   *     mask when properly set up)
   *   + VERTEX_OWNER: single-parameter macro returning the owner of a global
   *     vertex number
   *   + VERTEX_LOCAL: single-parameter macro returning the local offset of a
   *     global vertex number
   *   + VERTEX_TO_GLOBAL: single-parameter macro converting a local vertex
   *     offset to a global number
   *   + g.nlocalverts: number of vertices stored on the local rank
   *   + g.nglobalverts: total number of vertices in the graph
   *   + g.nlocaledges: number of graph edges stored locally
   *   + g.rowstarts, g.column: zero-based compressed sparse row data
   *     structure for the local part of the graph
   *
   * All macros documented above evaluate their arguments exactly once.
   *
   * The graph is stored using a 1-D, cyclic distribution: all edges incident
   * to vertex v are stored on rank (v % size) (aka VERTEX_OWNER(v)).  Edges
   * that are not self-loops are stored twice, once for each endpoint;
   * duplicates edges are kept.  The neighbors of vertex v can be obtained on
   * rank VERTEX_OWNER(v); they are stored in elements
   * {g.rowstarts[VERTEX_LOCAL(v)] ... g.rowstarts[VERTEX_LOCAL(v) + 1] - 1}
   * (inclusive) of g.column.
   *
   * Upon exit, your BFS must have filled in:
   *   + pred (an array of size g.nlocalverts):
   *     - The predecessor of vertex v in the BFS tree should go into
   *       pred[VERTEX_LOCAL(v)] on rank VERTEX_OWNER(v)
   *     - The predecessor of root is root
   *     - The predecessor of any unreachable vertex is -1
   *
   * The validator will check this for correctness. */
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
