/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include "onesided.h"
#include "common.h"
#include "../generator/make_graph.h"
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* These can be specific to this file, making them independent of whatever data
 * distribution the graph data structure and BFS are using. */

static inline int vertex_owner_for_find_roots(int64_t v) {
  return MOD_SIZE(v);
}

static inline size_t vertex_local_for_find_roots(int64_t v) {
  return DIV_SIZE(v);
}

static inline int64_t vertex_to_global_for_find_roots(int r, size_t l) {
  return MUL_SIZE((int64_t)l) + (int64_t)r;
}

/* Find num_bfs_roots random vertices, each of which has degree >= 1, with the
 * same answer produced on all nodes. */
void find_bfs_roots(int* num_bfs_roots, const int64_t nglobalverts, const tuple_graph* const tg, const uint64_t seed1, const uint64_t seed2, int64_t* const bfs_roots, int64_t* const max_used_vertex_ptr) {
  uint64_t nlocalverts = (nglobalverts / size) + (rank < nglobalverts % size);

  /* Create a bitmap (one byte per element) of vertices, setting those entries
   * that have nonzero degree. */
  unsigned char* restrict nonzero_degree = (unsigned char*)xmalloc(nlocalverts);
  memset(nonzero_degree, 0, nlocalverts);
  {
    unsigned char one = 1;
    int64_t max_bufsize = tuple_graph_max_bufsize(tg);
    scatter_constant* nonzero_degree_win = init_scatter_constant((void*)nonzero_degree, nlocalverts, sizeof(unsigned char), &one, 2 * size_min(HALF_CHUNKSIZE, max_bufsize), MPI_UNSIGNED_CHAR);
    ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize) {
      const packed_edge* restrict const edges = buf;
      ptrdiff_t ii;
      for (ii = 0; ii < max_bufsize; ii += HALF_CHUNKSIZE) {
        ptrdiff_t i;
        ptrdiff_t i_end = ptrdiff_min(ii + HALF_CHUNKSIZE, bufsize);
        begin_scatter_constant(nonzero_degree_win);
#pragma omp parallel for
        for (i = ii; i < i_end; ++i) {
          int64_t v0 = get_v0_from_edge(&edges[i]);
          int64_t v1 = get_v1_from_edge(&edges[i]);
          if (v0 == v1) continue;
          assert (v0 >= 0 && v0 < nglobalverts);
          assert (v1 >= 0 && v1 < nglobalverts);
          assert (vertex_owner_for_find_roots(v0) >= 0 && vertex_owner_for_find_roots(v0) < size);
          assert (vertex_owner_for_find_roots(v1) >= 0 && vertex_owner_for_find_roots(v1) < size);
          assert (vertex_owner_for_find_roots(v0) != rank || vertex_local_for_find_roots(v0) < nlocalverts);
          assert (vertex_owner_for_find_roots(v1) != rank || vertex_local_for_find_roots(v1) < nlocalverts);
          add_scatter_constant_request(nonzero_degree_win, vertex_owner_for_find_roots(v0), vertex_local_for_find_roots(v0), 2 * (i - ii) + 0);
          add_scatter_constant_request(nonzero_degree_win, vertex_owner_for_find_roots(v1), vertex_local_for_find_roots(v1), 2 * (i - ii) + 1);
        }
        end_scatter_constant(nonzero_degree_win);
      }
    } ITERATE_TUPLE_GRAPH_END;
    destroy_scatter_constant(nonzero_degree_win);
  }

  uint64_t counter = 0;
  int bfs_root_idx;
  for (bfs_root_idx = 0; bfs_root_idx < *num_bfs_roots; ++bfs_root_idx) {
    int64_t root;
    while (1) {
      double d[2];
      make_random_numbers(2, seed1, seed2, counter, d);
      root = (int64_t)((d[0] + d[1]) * nglobalverts) % nglobalverts;
      counter += 2;
      if (counter > 2 * nglobalverts) break;
      int is_duplicate = 0;
      int i;
      for (i = 0; i < bfs_root_idx; ++i) {
        if (root == bfs_roots[i]) {
          is_duplicate = 1;
          break;
        }
      }
      if (is_duplicate) continue; /* Everyone takes the same path here */
      int root_ok;
      if (vertex_owner_for_find_roots(root) == rank) {
        root_ok = nonzero_degree[vertex_local_for_find_roots(root)];
      }
      MPI_Bcast(&root_ok, 1, MPI_INT, vertex_owner_for_find_roots(root), MPI_COMM_WORLD);
      if (root_ok) break;
    }
    bfs_roots[bfs_root_idx] = root;
  }
  *num_bfs_roots = bfs_root_idx;

  /* Find maximum non-zero-degree vertex. */
  {
    size_t i;
    int64_t max_used_vertex = 0;
    for (i = nlocalverts; i > 0; --i) {
      if (nonzero_degree[i - 1]) {
        max_used_vertex = vertex_to_global_for_find_roots(rank, i - 1);
        break;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &max_used_vertex, 1, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);
    *max_used_vertex_ptr = max_used_vertex;
  }

  free(nonzero_degree);
}
