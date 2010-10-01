/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include "common.h"
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

/* Add your own BFS code into this file (or a copy of it). */
void run_mpi_bfs(const csr_graph* const g, int64_t root, int64_t* pred, int64_t* nvisited) {
  /* Predefined entities you can use in your BFS (from common.h):
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
   *   + g->nlocalverts: number of vertices stored on the local rank
   *   + g->nglobalverts: total number of vertices in the graph
   *   + g->nlocaledges: number of graph edges stored locally
   *   + g->rowstarts, g->column: zero-based compressed sparse row data
   *     structure for the local part of the graph
   *
   * All macros documented above evaluate their arguments exactly once.
   *
   * The graph is stored using a 1-D, cyclic distribution: all edges incident
   * to vertex v are stored on rank (v % size) (aka VERTEX_OWNER(v)).  Edges
   * that are not self-loops are stored twice, once for each endpoint;
   * duplicates edges are kept.  The neighbors of vertex v can be obtained on
   * rank VERTEX_OWNER(v); they are stored in elements
   * {g->rowstarts[VERTEX_LOCAL(v)] ... g->rowstarts[VERTEX_LOCAL(v) + 1] - 1}
   * (inclusive) of g->column.
   *
   * Upon exit, your BFS must have filled in:
   *   + pred (an array of size g->nlocalverts):
   *     - The predecessor of vertex v in the BFS tree should go into
   *       pred[VERTEX_LOCAL(v)] on rank VERTEX_OWNER(v)
   *     - The predecessor of root is root
   *     - The predecessor of any unreachable vertex is -1
   *   + *nvisited
   *     - The number of vertices, including the root, visited during the BFS
   *     - This number should be summed over all ranks and should be consistent
   *
   * The validator will check both of these for correctness. */
}
