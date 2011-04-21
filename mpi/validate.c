/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include "onesided.h"
#include "common.h"
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

int validate_bfs_result(const tuple_graph* const tg, const int64_t nglobalverts, const size_t nlocalverts, const int64_t root, const int64_t* const pred, int64_t* const edge_visit_count_ptr) {
  int validation_passed = 1;
  int root_owner;
  size_t root_local;
  get_vertex_distribution_for_pred(1, &root, &root_owner, &root_local);
  int root_is_mine = (root_owner == rank);

  /* Get maximum values so loop counts are consistent across ranks. */
  uint64_t maxlocalverts_ui = nlocalverts;
  MPI_Allreduce(MPI_IN_PLACE, &maxlocalverts_ui, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  size_t maxlocalverts = (size_t)maxlocalverts_ui;

  ptrdiff_t max_bufsize = tuple_graph_max_bufsize(tg);
  ptrdiff_t edge_chunk_size = ptrdiff_min(HALF_CHUNKSIZE, max_bufsize);

  /* Check that root is its own parent. */
  if (root_is_mine) {
    if (pred[root_local] != root) {
      fprintf(stderr, "%d: Validation error: parent of root vertex %" PRId64 " is %" PRId64 ", not the root itself.\n", rank, root, pred[root_local]);
      validation_passed = 0;
    }
  }

  /* Check that nothing else is its own parent, and check for in-range
   * values. */
  int any_range_errors = 0;
  {
    int* restrict pred_owner = (int*)xmalloc(size_min(CHUNKSIZE, nlocalverts) * sizeof(int));
    size_t* restrict pred_local = (size_t*)xmalloc(size_min(CHUNKSIZE, nlocalverts) * sizeof(size_t));
    size_t ii;
    for (ii = 0; ii < nlocalverts; ii += CHUNKSIZE) {
      ptrdiff_t i_end = ptrdiff_min(ii + CHUNKSIZE, nlocalverts);
      get_vertex_distribution_for_pred(i_end - ii, &pred[ii], pred_owner, pred_local);
      ptrdiff_t i;
#pragma omp parallel for reduction(&&:validation_passed) reduction(||:any_range_errors)
      for (i = ii; i < i_end; ++i) {
        if ((!root_is_mine || i != root_local) && pred[i] != -1 && pred_owner[i - ii] == rank && pred_local[i - ii] == i) {
          fprintf(stderr, "%d: Validation error: parent of non-root vertex %" PRId64 " is itself.\n", rank, vertex_to_global_for_pred(rank, i));
          validation_passed = 0;
        }
        if (pred[i] < -1 || pred[i] >= nglobalverts) {
          fprintf(stderr, "%d: Validation error: parent of vertex %" PRId64 " is out-of-range value %" PRId64 ".\n", rank, vertex_to_global_for_pred(rank, i), pred[i]);
          validation_passed = 0;
          any_range_errors = 1;
        }
      }
    }
    free(pred_owner);
    free(pred_local);
  }
  MPI_Allreduce(MPI_IN_PLACE, &any_range_errors, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

  if (!any_range_errors) { /* Other parts of validation assume in-range values */

    /* Create a vertex depth map to use for later validation. */
    uint16_t* restrict depth = (uint16_t*)xMPI_Alloc_mem(nlocalverts * sizeof(uint16_t));
    {
      {
        ptrdiff_t i;
#pragma omp parallel for
        for (i = 0; i < nlocalverts; ++i) depth[i] = UINT16_MAX;
        if (root_is_mine) depth[root_local] = 0;
      }
      uint16_t* restrict pred_depth = (uint16_t*)xMPI_Alloc_mem(size_min(CHUNKSIZE, nlocalverts) * sizeof(uint16_t)); /* Depth of predecessor vertex for each local vertex */
      gather* depth_win = init_gather((void*)depth, nlocalverts, sizeof(uint16_t), pred_depth, size_min(CHUNKSIZE, nlocalverts), MPI_UINT16_T);
      int* restrict pred_owner = (int*)xmalloc(size_min(CHUNKSIZE, nlocalverts) * sizeof(int));
      size_t* restrict pred_local = (size_t*)xmalloc(size_min(CHUNKSIZE, nlocalverts) * sizeof(size_t));
      {
        /* Iteratively update depth[v] = min(depth[v], depth[pred[v]] + 1) [saturating at UINT16_MAX] until no changes. */
        while (1) {
#if 0
          fprintf(stderr, "Rank %d: depth map iteration\n", rank);
#endif
          int any_changes = 0;
          size_t ii;
          for (ii = 0; ii < maxlocalverts; ii += CHUNKSIZE) {
            ptrdiff_t i_end = ptrdiff_min(ii + CHUNKSIZE, nlocalverts);
            begin_gather(depth_win);
            get_vertex_distribution_for_pred(i_end - ii, &pred[ii], pred_owner, pred_local);
            ptrdiff_t i;
#pragma omp parallel for
            for (i = ii; i < i_end; ++i) {
              if (pred[i] != -1) {
                add_gather_request(depth_win, i - ii, pred_owner[i - ii], pred_local[i - ii], i - ii);
              } else {
                pred_depth[i - ii] = UINT16_MAX;
              }
            }
            end_gather(depth_win);
#pragma omp parallel for reduction(&&:validation_passed) reduction(||:any_changes)
            for (i = ii; i < i_end; ++i) {
              if (rank == root_owner && i == root_local) continue;
              if (pred_depth[i - ii] != UINT16_MAX) {
                if (depth[i] != UINT16_MAX && depth[i] != pred_depth[i - ii] + 1) {
                  fprintf(stderr, "%d: Validation error: BFS predecessors do not form a tree; see vertices %" PRId64 " (depth %" PRIu16 ") and %" PRId64 " (depth %" PRIu16 ").\n", rank, vertex_to_global_for_pred(rank, i), depth[i], pred[i - ii], pred_depth[i - ii]);
                  validation_passed = 0;
                } else if (depth[i] == pred_depth[i - ii] + 1) {
                  /* Nothing to do */
                } else {
                  depth[i] = pred_depth[i - ii] + 1;
                  any_changes = 1;
                }
              }
            }
          }
          MPI_Allreduce(MPI_IN_PLACE, &any_changes, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
          if (!any_changes) break;
        }
      }
      destroy_gather(depth_win);
      MPI_Free_mem(pred_depth);
      free(pred_owner);
      free(pred_local);
    }

    {
      /* Check that all edges connect vertices whose depths differ by at most
       * one, and check that there is an edge from each vertex to its claimed
       * predecessor.  Also, count visited edges (including duplicates and
       * self-loops).  */
      unsigned char* restrict pred_valid = (unsigned char*)xMPI_Alloc_mem(nlocalverts * sizeof(unsigned char));
      memset(pred_valid, 0, nlocalverts * sizeof(unsigned char));
      int64_t* restrict edge_endpoint = (int64_t*)xmalloc(2 * edge_chunk_size * sizeof(int64_t));
      int* restrict edge_owner = (int*)xmalloc(2 * edge_chunk_size * sizeof(int));
      size_t* restrict edge_local = (size_t*)xmalloc(2 * edge_chunk_size * sizeof(size_t));
      int64_t* restrict edge_preds = (int64_t*)xMPI_Alloc_mem(2 * edge_chunk_size * sizeof(int64_t));
      gather* pred_win = init_gather((void*)pred, nlocalverts, sizeof(int64_t), edge_preds, 2 * edge_chunk_size, MPI_INT64_T);
      unsigned char one = 1;
      scatter_constant* pred_valid_win = init_scatter_constant((void*)pred_valid, nlocalverts, sizeof(unsigned char), &one, 2 * edge_chunk_size, MPI_UNSIGNED_CHAR);
      int64_t edge_visit_count = 0;
      uint16_t* restrict edge_depths = (uint16_t*)xMPI_Alloc_mem(2 * edge_chunk_size * sizeof(uint16_t));
      gather* depth_win = init_gather((void*)depth, nlocalverts, sizeof(uint16_t), edge_depths, 2 * edge_chunk_size, MPI_UINT16_T);
      ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize) {
        ptrdiff_t ii;
        for (ii = 0; ii < max_bufsize; ii += HALF_CHUNKSIZE) {
          ptrdiff_t i_end = ptrdiff_min(ii + HALF_CHUNKSIZE, bufsize);
          ptrdiff_t i;
#pragma omp parallel for
          for (i = ii; i < i_end; ++i) {
            int64_t v0 = get_v0_from_edge(&buf[i]);
            int64_t v1 = get_v1_from_edge(&buf[i]);
            edge_endpoint[(i - ii) * 2 + 0] = v0;
            edge_endpoint[(i - ii) * 2 + 1] = v1;
          }
          get_vertex_distribution_for_pred(2 * (i_end - ii), edge_endpoint, edge_owner, edge_local);
          begin_gather(depth_win);
          begin_gather(pred_win);
#pragma omp parallel for
          for (i = ii; i < i_end; ++i) {
            add_gather_request(depth_win, (i - ii) * 2 + 0, edge_owner[(i - ii) * 2 + 0], edge_local[(i - ii) * 2 + 0], (i - ii) * 2 + 0);
            add_gather_request(depth_win, (i - ii) * 2 + 1, edge_owner[(i - ii) * 2 + 1], edge_local[(i - ii) * 2 + 1], (i - ii) * 2 + 1);
            add_gather_request(pred_win, (i - ii) * 2 + 0, edge_owner[(i - ii) * 2 + 0], edge_local[(i - ii) * 2 + 0], (i - ii) * 2 + 0);
            add_gather_request(pred_win, (i - ii) * 2 + 1, edge_owner[(i - ii) * 2 + 1], edge_local[(i - ii) * 2 + 1], (i - ii) * 2 + 1);
          }
          end_gather(pred_win);
          end_gather(depth_win);
          begin_scatter_constant(pred_valid_win);
#pragma omp parallel for reduction(&&:validation_passed) reduction(+:edge_visit_count)
          for (i = ii; i < i_end; ++i) {
            int64_t src = get_v0_from_edge(&buf[i]);
            int64_t tgt = get_v1_from_edge(&buf[i]);
            uint16_t src_depth = edge_depths[(i - ii) * 2 + 0];
            uint16_t tgt_depth = edge_depths[(i - ii) * 2 + 1];
            if (src_depth != UINT16_MAX && tgt_depth == UINT16_MAX) {
              fprintf(stderr, "%d: Validation error: edge connects vertex %" PRId64 " in the BFS tree (depth %" PRIu16 ") to vertex %" PRId64 " outside the tree.\n", rank, src, src_depth, tgt);
              validation_passed = 0;
            } else if (src_depth == UINT16_MAX && tgt_depth != UINT16_MAX) {
              fprintf(stderr, "%d: Validation error: edge connects vertex %" PRId64 " in the BFS tree (depth %" PRIu16 ") to vertex %" PRId64 " outside the tree.\n", rank, tgt, tgt_depth, src);
              validation_passed = 0;
            } else if (src_depth - tgt_depth < -1 ||
                       src_depth - tgt_depth > 1) {
              fprintf(stderr, "%d: Validation error: depths of edge endpoints %" PRId64 " (depth %" PRIu16 ") and %" PRId64 " (depth %" PRIu16 ") are too far apart (abs. val. > 1).\n", rank, src, src_depth, tgt, tgt_depth);
              validation_passed = 0;
            } else if (src_depth != UINT16_MAX) {
              ++edge_visit_count;
            }
            if (edge_preds[(i - ii) * 2 + 0] == tgt) {
              add_scatter_constant_request(pred_valid_win, edge_owner[(i - ii) * 2 + 0], edge_local[(i - ii) * 2 + 0], (i - ii) * 2 + 0);
            }
            if (edge_preds[(i - ii) * 2 + 1] == src) {
              add_scatter_constant_request(pred_valid_win, edge_owner[(i - ii) * 2 + 1], edge_local[(i - ii) * 2 + 1], (i - ii) * 2 + 1);
            }
          }
          end_scatter_constant(pred_valid_win);
        }
      } ITERATE_TUPLE_GRAPH_END;
      destroy_gather(depth_win);
      MPI_Free_mem(edge_depths);
      destroy_gather(pred_win);
      MPI_Free_mem(edge_preds);
      free(edge_owner);
      free(edge_local);
      free(edge_endpoint);
      destroy_scatter_constant(pred_valid_win);
      ptrdiff_t i;
#pragma omp parallel for reduction(&&:validation_passed)
      for (i = 0; i < nlocalverts; ++i) {
        int64_t p = pred[i];
        if (p == -1) continue;
        int found_pred_edge = pred_valid[i];
        if (root_owner == rank && root_local == i) found_pred_edge = 1; /* Root vertex */
        if (!found_pred_edge) {
          int64_t v = vertex_to_global_for_pred(rank, i);
          fprintf(stderr, "%d: Validation error: no graph edge from vertex %" PRId64 " to its parent %" PRId64 ".\n", rank, v, pred[i]);
          validation_passed = 0;
        }
      }
      MPI_Free_mem(pred_valid);

      MPI_Allreduce(MPI_IN_PLACE, &edge_visit_count, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
      *edge_visit_count_ptr = edge_visit_count;
    }
    MPI_Free_mem(depth);

  } /* End of part skipped by range errors */
  
  /* Collect the global validation result. */
  MPI_Allreduce(MPI_IN_PLACE, &validation_passed, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  return validation_passed;
}
