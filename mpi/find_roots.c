/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include "common.h"
#include "../generator/make_graph.h"
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

/* Find num_bfs_roots random vertices, each of which has degree >= 1, with the
 * same answer produced on all nodes. */
void find_bfs_roots(int* num_bfs_roots, const csr_graph* const g, const uint64_t seed1, const uint64_t seed2, int64_t* const bfs_roots) {
  /* This implementation is slow, but there aren't enough roots being
   * generated for that to be a big issue. */
  uint64_t counter = 0;
  int bfs_root_idx;
  for (bfs_root_idx = 0; bfs_root_idx < *num_bfs_roots; ++bfs_root_idx) {
    int64_t root;
    while (1) {
      double d[2];
      make_random_numbers(2, seed1, seed2, counter, d);
      root = (int64_t)((d[0] + d[1]) * g->nglobalverts) % g->nglobalverts;
      counter += 2;
      if (counter > 2*g->nglobalverts) break;
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
      if (VERTEX_OWNER(root) == rank) {
        root_ok = 0;
        size_t ei, ei_end = g->rowstarts[VERTEX_LOCAL(root) + 1];
        for (ei = g->rowstarts[VERTEX_LOCAL(root)]; ei < ei_end; ++ei) {
          if (g->column[ei] != root) {
            root_ok = 1;
            break;
          }
        }
      }
      MPI_Bcast(&root_ok, 1, MPI_INT, VERTEX_OWNER(root), MPI_COMM_WORLD);
      if (root_ok) break;
    }
    bfs_roots[bfs_root_idx] = root;
  }
  *num_bfs_roots = bfs_root_idx;
}
