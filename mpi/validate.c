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

/* This code assumes signed shifts are arithmetic, which they are on
 * practically all modern systems but is not guaranteed by C. */

static inline int64_t get_pred_from_pred_entry(int64_t val) {
  return (val << 16) >> 16;
}

static inline uint16_t get_depth_from_pred_entry(int64_t val) {
  return (val >> 48) & 0xFFFF;
}

static inline void write_pred_entry_depth(int64_t* loc, uint16_t depth) {
  *loc = (*loc & INT64_C(0xFFFFFFFFFFFF)) | ((int64_t)(depth & 0xFFFF) << 48);
}

/* Returns true if all values are in range. */
static int check_value_ranges(const int64_t nglobalverts, const size_t nlocalverts, const int64_t* const pred) {
  int any_range_errors = 0;
  {
    size_t ii;
    for (ii = 0; ii < nlocalverts; ii += CHUNKSIZE) {
      ptrdiff_t i_start = ii;
      ptrdiff_t i_end = ptrdiff_min(ii + CHUNKSIZE, nlocalverts);
      ptrdiff_t i;
      assert (i_start >= 0 && i_start <= (ptrdiff_t)nlocalverts);
      assert (i_end >= 0 && i_end <= (ptrdiff_t)nlocalverts);
#pragma omp parallel for reduction(||:any_range_errors)
      for (i = i_start; i < i_end; ++i) {
        int64_t p = get_pred_from_pred_entry(pred[i]);
        if (p < -1 || p >= nglobalverts) {
          fprintf(stderr, "%d: Validation error: parent of vertex %" PRId64 " is out-of-range value %" PRId64 ".\n", rank, vertex_to_global_for_pred(rank, i), p);
          any_range_errors = 1;
        }
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &any_range_errors, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
  return !any_range_errors;
}

/* Use the predecessors in the given map to write the BFS levels to the high 16
 * bits of each element in pred; this also catches some problems in pred
 * itself.  Returns true if the predecessor map is valid. */
static int build_bfs_depth_map(const int64_t nglobalverts, const size_t nlocalverts, const size_t maxlocalverts, const int64_t root, int64_t* const pred) {
  (void)nglobalverts;
  int validation_passed = 1;
  int root_owner;
  size_t root_local;
  get_vertex_distribution_for_pred(1, &root, &root_owner, &root_local);
  int root_is_mine = (root_owner == rank);
  if (root_is_mine) assert (root_local < nlocalverts);

  {
    ptrdiff_t i;
#pragma omp parallel for
    for (i = 0; i < (ptrdiff_t)nlocalverts; ++i) write_pred_entry_depth(&pred[i], UINT16_MAX);
    if (root_is_mine) write_pred_entry_depth(&pred[root_local], 0);
  }
  int64_t* restrict pred_pred = (int64_t*)xMPI_Alloc_mem(size_min(CHUNKSIZE, nlocalverts) * sizeof(int64_t)); /* Predecessor info of predecessor vertex for each local vertex */
  gather* pred_win = init_gather((void*)pred, nlocalverts, sizeof(int64_t), pred_pred, size_min(CHUNKSIZE, nlocalverts), size_min(CHUNKSIZE, nlocalverts), MPI_INT64_T);
  int64_t* restrict pred_vtx = (int64_t*)xmalloc(size_min(CHUNKSIZE, nlocalverts) * sizeof(int64_t)); /* Vertex (not depth) part of pred map */
  int* restrict pred_owner = (int*)xmalloc(size_min(CHUNKSIZE, nlocalverts) * sizeof(int));
  size_t* restrict pred_local = (size_t*)xmalloc(size_min(CHUNKSIZE, nlocalverts) * sizeof(size_t));
  int iter_number = 0;
  {
    /* Iteratively update depth[v] = min(depth[v], depth[pred[v]] + 1) [saturating at UINT16_MAX] until no changes. */
    while (1) {
      ++iter_number;
      int any_changes = 0;
      ptrdiff_t ii;
      for (ii = 0; ii < (ptrdiff_t)maxlocalverts; ii += CHUNKSIZE) {
        ptrdiff_t i_start = ptrdiff_min(ii, nlocalverts);
        ptrdiff_t i_end = ptrdiff_min(ii + CHUNKSIZE, nlocalverts);
        begin_gather(pred_win);
        ptrdiff_t i;
        assert (i_start >= 0 && i_start <= (ptrdiff_t)nlocalverts);
        assert (i_end >= 0 && i_end <= (ptrdiff_t)nlocalverts);
#pragma omp parallel for
        for (i = i_start; i < i_end; ++i) {
          pred_vtx[i - i_start] = get_pred_from_pred_entry(pred[i]);
        }
        get_vertex_distribution_for_pred(i_end - i_start, pred_vtx, pred_owner, pred_local);
#pragma omp parallel for
        for (i = i_start; i < i_end; ++i) {
          if (pred[i] != -1) {
            add_gather_request(pred_win, i - i_start, pred_owner[i - i_start], pred_local[i - i_start], i - i_start);
          } else {
            pred_pred[i - i_start] = -1;
          }
        }
        end_gather(pred_win);
#pragma omp parallel for reduction(&&:validation_passed) reduction(||:any_changes)
        for (i = i_start; i < i_end; ++i) {
          if (rank == root_owner && (size_t)i == root_local) continue;
          if (get_depth_from_pred_entry(pred_pred[i - i_start]) != UINT16_MAX) {
            if (get_depth_from_pred_entry(pred[i]) != UINT16_MAX && get_depth_from_pred_entry(pred[i]) != get_depth_from_pred_entry(pred_pred[i - i_start]) + 1) {
              fprintf(stderr, "%d: Validation error: BFS predecessors do not form a tree; see vertices %" PRId64 " (depth %" PRIu16 ") and %" PRId64 " (depth %" PRIu16 ").\n", rank, vertex_to_global_for_pred(rank, i), get_depth_from_pred_entry(pred[i]), get_pred_from_pred_entry(pred[i]), get_depth_from_pred_entry(pred_pred[i - i_start]));
              validation_passed = 0;
            } else if (get_depth_from_pred_entry(pred[i]) == get_depth_from_pred_entry(pred_pred[i - i_start]) + 1) {
              /* Nothing to do */
            } else {
              write_pred_entry_depth(&pred[i], get_depth_from_pred_entry(pred_pred[i - i_start]) + 1);
              any_changes = 1;
            }
          }
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &any_changes, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
      if (!any_changes) break;
    }
  }
  destroy_gather(pred_win);
  MPI_Free_mem(pred_pred);
  free(pred_owner);
  free(pred_local);
  free(pred_vtx);
  return validation_passed;
}

/* Check the BFS levels in pred against the predecessors given there.  Returns
 * true if the maps are valid. */
static int check_bfs_depth_map_using_predecessors(const tuple_graph* const tg, const int64_t nglobalverts, const size_t nlocalverts, const size_t maxlocalverts, const int64_t root, const int64_t* const pred) {
  (void)nglobalverts; /* Avoid warning */
  assert (tg->edgememory_size >= 0 && tg->max_edgememory_size >= tg->edgememory_size && tg->max_edgememory_size <= tg->nglobaledges);
  assert (root >= 0 && root < nglobalverts);
  assert (nglobalverts >= 0);
  assert (pred);

  int validation_passed = 1;
  int root_owner;
  size_t root_local;
  get_vertex_distribution_for_pred(1, &root, &root_owner, &root_local);
  int root_is_mine = (root_owner == rank);
  if (root_is_mine) assert (root_local < nlocalverts);

  {
    ptrdiff_t i;
    if (root_is_mine && get_depth_from_pred_entry(pred[root_local]) != 0) {
      fprintf(stderr, "%d: Validation error: depth of root vertex %" PRId64 " is %" PRIu16 ", not 0.\n", rank, root, get_depth_from_pred_entry(pred[root_local]));
      validation_passed = 0;
    }
#pragma omp parallel for reduction(&&:validation_passed)
    for (i = 0; i < (ptrdiff_t)nlocalverts; ++i) {
      if (get_pred_from_pred_entry(pred[i]) == -1 &&
          get_depth_from_pred_entry(pred[i]) != UINT16_MAX) {
        fprintf(stderr, "%d: Validation error: depth of vertex %" PRId64 " with no predecessor is %" PRIu16 ", not UINT16_MAX.\n", rank, vertex_to_global_for_pred(rank, i), get_depth_from_pred_entry(pred[i]));
        validation_passed = 0;
      } else if (get_pred_from_pred_entry(pred[i]) != -1 &&
                 get_depth_from_pred_entry(pred[i]) == UINT16_MAX) {
        fprintf(stderr, "%d: Validation error: predecessor of claimed unreachable vertex %" PRId64 " is %" PRId64 ", not -1.\n", rank, vertex_to_global_for_pred(rank, i), get_pred_from_pred_entry(pred[i]));
        validation_passed = 0;
      }
    }
  }
  int64_t* restrict pred_pred = (int64_t*)xMPI_Alloc_mem(size_min(CHUNKSIZE, nlocalverts) * sizeof(int64_t)); /* Predecessor info of predecessor vertex for each local vertex */
  gather* pred_win = init_gather((void*)pred, nlocalverts, sizeof(int64_t), pred_pred, size_min(CHUNKSIZE, nlocalverts), size_min(CHUNKSIZE, nlocalverts), MPI_INT64_T);
  int64_t* restrict pred_vtx = (int64_t*)xmalloc(size_min(CHUNKSIZE, nlocalverts) * sizeof(int64_t)); /* Vertex (not depth) part of pred map */
  int* restrict pred_owner = (int*)xmalloc(size_min(CHUNKSIZE, nlocalverts) * sizeof(int));
  size_t* restrict pred_local = (size_t*)xmalloc(size_min(CHUNKSIZE, nlocalverts) * sizeof(size_t));
  size_t ii;
  for (ii = 0; ii < maxlocalverts; ii += CHUNKSIZE) {
    ptrdiff_t i_start = ptrdiff_min(ii, nlocalverts);
    ptrdiff_t i_end = ptrdiff_min(ii + CHUNKSIZE, nlocalverts);
    begin_gather(pred_win);
    ptrdiff_t i;
    assert (i_start >= 0 && i_start <= (ptrdiff_t)nlocalverts);
    assert (i_end >= 0 && i_end <= (ptrdiff_t)nlocalverts);
    assert (i_end >= i_start);
    assert (i_end - i_start >= 0 && i_end - i_start <= (ptrdiff_t)size_min(CHUNKSIZE, nlocalverts));
#pragma omp parallel for
    for (i = i_start; i < i_end; ++i) {
      pred_vtx[i - i_start] = get_pred_from_pred_entry(pred[i]);
    }
    get_vertex_distribution_for_pred(i_end - i_start, pred_vtx, pred_owner, pred_local);
#pragma omp parallel for
    for (i = i_start; i < i_end; ++i) {
      if (pred[i] != -1) {
        add_gather_request(pred_win, i - i_start, pred_owner[i - i_start], pred_local[i - i_start], i - i_start);
      } else {
        pred_pred[i - i_start] = -1;
      }
    }
    end_gather(pred_win);
#pragma omp parallel for reduction(&&:validation_passed)
    for (i = i_start; i < i_end; ++i) {
      if (rank == root_owner && (size_t)i == root_local) continue;
      if (get_pred_from_pred_entry(pred[i]) == -1) continue; /* Already checked */
      if (get_depth_from_pred_entry(pred_pred[i - i_start]) == UINT16_MAX) {
        fprintf(stderr, "%d: Validation error: predecessor %" PRId64 " of vertex %" PRId64 " (depth %" PRIu16 ") is marked as unreachable.\n", rank, get_pred_from_pred_entry(pred[i]), vertex_to_global_for_pred(rank, i), get_depth_from_pred_entry(pred[i]));
        validation_passed = 0;
      }
      if (get_depth_from_pred_entry(pred[i]) != get_depth_from_pred_entry(pred_pred[i - i_start]) + 1) {
        fprintf(stderr, "%d: Validation error: BFS predecessors do not form a tree; see vertices %" PRId64 " (depth %" PRIu16 ") and %" PRId64 " (depth %" PRIu16 ").\n", rank, vertex_to_global_for_pred(rank, i), get_depth_from_pred_entry(pred[i]), get_pred_from_pred_entry(pred[i]), get_depth_from_pred_entry(pred_pred[i - i_start]));
        validation_passed = 0;
      }
    }
  }
  destroy_gather(pred_win);
  MPI_Free_mem(pred_pred);
  free(pred_owner);
  free(pred_local);
  free(pred_vtx);
  return validation_passed;
}
/* Returns true if result is valid.  Also, updates high 16 bits of each element
 * of pred to contain the BFS level number (or -1 if not visited) of each
 * vertex; this is based on the predecessor map if the user didn't provide it.
 * */
int validate_bfs_result(const tuple_graph* const tg, const int64_t nglobalverts, const size_t nlocalverts, const int64_t root, int64_t* const pred, int64_t* const edge_visit_count_ptr) {

  assert (tg->edgememory_size >= 0 && tg->max_edgememory_size >= tg->edgememory_size && tg->max_edgememory_size <= tg->nglobaledges);
  assert (pred);
  *edge_visit_count_ptr = 0; /* Ensure it is a valid pointer */
  int ranges_ok = check_value_ranges(nglobalverts, nlocalverts, pred);
  if (root < 0 || root >= nglobalverts) {
    fprintf(stderr, "%d: Validation error: root vertex %" PRId64 " is invalid.\n", rank, root);
    ranges_ok = 0;
  }
  if (!ranges_ok) return 0; /* Fail */

  assert (tg->edgememory_size >= 0 && tg->max_edgememory_size >= tg->edgememory_size && tg->max_edgememory_size <= tg->nglobaledges);
  assert (pred);

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

  assert (tg->edgememory_size >= 0 && tg->max_edgememory_size >= tg->edgememory_size && tg->max_edgememory_size <= tg->nglobaledges);
  assert (pred);

  /* Check that root is its own parent. */
  if (root_is_mine) {
    assert (root_local < nlocalverts);
    if (get_pred_from_pred_entry(pred[root_local]) != root) {
      fprintf(stderr, "%d: Validation error: parent of root vertex %" PRId64 " is %" PRId64 ", not the root itself.\n", rank, root, get_pred_from_pred_entry(pred[root_local]));
      validation_passed = 0;
    }
  }

  assert (tg->edgememory_size >= 0 && tg->max_edgememory_size >= tg->edgememory_size && tg->max_edgememory_size <= tg->nglobaledges);
  assert (pred);

  /* Check that nothing else is its own parent. */
  {
    int* restrict pred_owner = (int*)xmalloc(size_min(CHUNKSIZE, nlocalverts) * sizeof(int));
    size_t* restrict pred_local = (size_t*)xmalloc(size_min(CHUNKSIZE, nlocalverts) * sizeof(size_t));
    int64_t* restrict pred_vtx = (int64_t*)xmalloc(size_min(CHUNKSIZE, nlocalverts) * sizeof(int64_t)); /* Vertex (not depth) part of pred map */
    ptrdiff_t ii;
    for (ii = 0; ii < (ptrdiff_t)nlocalverts; ii += CHUNKSIZE) {
      ptrdiff_t i_start = ii;
      ptrdiff_t i_end = ptrdiff_min(ii + CHUNKSIZE, nlocalverts);
      ptrdiff_t i;
      assert (i_start >= 0 && i_start <= (ptrdiff_t)nlocalverts);
      assert (i_end >= 0 && i_end <= (ptrdiff_t)nlocalverts);
#pragma omp parallel for
      for (i = i_start; i < i_end; ++i) {
        pred_vtx[i - i_start] = get_pred_from_pred_entry(pred[i]);
      }
      get_vertex_distribution_for_pred(i_end - i_start, pred_vtx, pred_owner, pred_local);
#pragma omp parallel for reduction(&&:validation_passed)
      for (i = i_start; i < i_end; ++i) {
        if ((!root_is_mine || (size_t)i != root_local) &&
            get_pred_from_pred_entry(pred[i]) != -1 &&
            pred_owner[i - i_start] == rank &&
            pred_local[i - i_start] == (size_t)i) {
          fprintf(stderr, "%d: Validation error: parent of non-root vertex %" PRId64 " is itself.\n", rank, vertex_to_global_for_pred(rank, i));
          validation_passed = 0;
        }
      }
    }
    free(pred_owner);
    free(pred_local);
    free(pred_vtx);
  }

  assert (tg->edgememory_size >= 0 && tg->max_edgememory_size >= tg->edgememory_size && tg->max_edgememory_size <= tg->nglobaledges);
  assert (pred);

  if (bfs_writes_depth_map()) {
    int check_ok = check_bfs_depth_map_using_predecessors(tg, nglobalverts, nlocalverts, maxlocalverts, root, pred);
    if (!check_ok) validation_passed = 0;
  } else {
    /* Create a vertex depth map to use for later validation. */
    int pred_ok = build_bfs_depth_map(nglobalverts, nlocalverts, maxlocalverts, root, pred);
    if (!pred_ok) validation_passed = 0;
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
    gather* pred_win = init_gather((void*)pred, nlocalverts, sizeof(int64_t), edge_preds, 2 * edge_chunk_size, 2 * edge_chunk_size, MPI_INT64_T);
    unsigned char one = 1;
    scatter_constant* pred_valid_win = init_scatter_constant((void*)pred_valid, nlocalverts, sizeof(unsigned char), &one, 2 * edge_chunk_size, MPI_UNSIGNED_CHAR);
    int64_t edge_visit_count = 0;
    ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize) {
      ptrdiff_t ii;
      for (ii = 0; ii < max_bufsize; ii += HALF_CHUNKSIZE) {
        ptrdiff_t i_start = ptrdiff_min(ii, bufsize);
        ptrdiff_t i_end = ptrdiff_min(ii + HALF_CHUNKSIZE, bufsize);
        assert (i_end - i_start <= edge_chunk_size);
        ptrdiff_t i;
#pragma omp parallel for
        for (i = i_start; i < i_end; ++i) {
          int64_t v0 = get_v0_from_edge(&buf[i]);
          int64_t v1 = get_v1_from_edge(&buf[i]);
          edge_endpoint[(i - i_start) * 2 + 0] = v0;
          edge_endpoint[(i - i_start) * 2 + 1] = v1;
        }
        get_vertex_distribution_for_pred(2 * (i_end - i_start), edge_endpoint, edge_owner, edge_local);
        begin_gather(pred_win);
#pragma omp parallel for
        for (i = i_start; i < i_end; ++i) {
          add_gather_request(pred_win, (i - i_start) * 2 + 0, edge_owner[(i - i_start) * 2 + 0], edge_local[(i - i_start) * 2 + 0], (i - i_start) * 2 + 0);
          add_gather_request(pred_win, (i - i_start) * 2 + 1, edge_owner[(i - i_start) * 2 + 1], edge_local[(i - i_start) * 2 + 1], (i - i_start) * 2 + 1);
        }
        end_gather(pred_win);
        begin_scatter_constant(pred_valid_win);
#pragma omp parallel for reduction(&&:validation_passed) reduction(+:edge_visit_count)
        for (i = i_start; i < i_end; ++i) {
          int64_t src = get_v0_from_edge(&buf[i]);
          int64_t tgt = get_v1_from_edge(&buf[i]);
          uint16_t src_depth = get_depth_from_pred_entry(edge_preds[(i - i_start) * 2 + 0]);
          uint16_t tgt_depth = get_depth_from_pred_entry(edge_preds[(i - i_start) * 2 + 1]);
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
          if (get_pred_from_pred_entry(edge_preds[(i - i_start) * 2 + 0]) == tgt) {
            add_scatter_constant_request(pred_valid_win, edge_owner[(i - i_start) * 2 + 0], edge_local[(i - i_start) * 2 + 0], (i - i_start) * 2 + 0);
          }
          if (get_pred_from_pred_entry(edge_preds[(i - i_start) * 2 + 1]) == src) {
            add_scatter_constant_request(pred_valid_win, edge_owner[(i - i_start) * 2 + 1], edge_local[(i - i_start) * 2 + 1], (i - i_start) * 2 + 1);
          }
        }
        end_scatter_constant(pred_valid_win);
      }
    } ITERATE_TUPLE_GRAPH_END;
    destroy_gather(pred_win);
    MPI_Free_mem(edge_preds);
    free(edge_owner);
    free(edge_local);
    free(edge_endpoint);
    destroy_scatter_constant(pred_valid_win);
    ptrdiff_t i;
#pragma omp parallel for reduction(&&:validation_passed)
    for (i = 0; i < (ptrdiff_t)nlocalverts; ++i) {
      int64_t p = get_pred_from_pred_entry(pred[i]);
      if (p == -1) continue;
      int found_pred_edge = pred_valid[i];
      if (root_owner == rank && root_local == (size_t)i) found_pred_edge = 1; /* Root vertex */
      if (!found_pred_edge) {
        int64_t v = vertex_to_global_for_pred(rank, i);
        fprintf(stderr, "%d: Validation error: no graph edge from vertex %" PRId64 " to its parent %" PRId64 ".\n", rank, v, get_pred_from_pred_entry(pred[i]));
        validation_passed = 0;
      }
    }
    MPI_Free_mem(pred_valid);

    MPI_Allreduce(MPI_IN_PLACE, &edge_visit_count, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    *edge_visit_count_ptr = edge_visit_count;
  }

  /* Collect the global validation result. */
  MPI_Allreduce(MPI_IN_PLACE, &validation_passed, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  return validation_passed;
}
