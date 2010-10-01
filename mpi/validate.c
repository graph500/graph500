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
#include <stdio.h>
#include <assert.h>

int validate_bfs_result(const csr_graph* const g, const int64_t root, const int64_t* const pred, const int64_t nvisited) {
  int validation_passed = 1;
  int root_is_mine = (VERTEX_OWNER(root) == rank);

  const size_t nlocalverts = g->nlocalverts;
  const size_t nlocaledges = g->nlocaledges;
  const int64_t nglobalverts = g->nglobalverts;

  /* Check that root is its own parent. */
  if (root_is_mine) {
    if (pred[VERTEX_LOCAL(root)] != root) {
      fprintf(stderr, "%d: Validation error: parent of root vertex %" PRId64 " is %" PRId64 ", not the root itself.\n", rank, root, pred[VERTEX_LOCAL(root)]);
      validation_passed = 0;
    }
  }

  /* Check that nothing else is its own parent, and check for in-range
   * values. */
  int any_range_errors = 0;
  size_t i;
  for (i = 0; i < nlocalverts; ++i) {
    int64_t v = VERTEX_TO_GLOBAL(i);
    assert (VERTEX_OWNER(v) == rank);
    assert (VERTEX_LOCAL(v) == i);
    if (v != root && pred[i] == v) {
      fprintf(stderr, "%d: Validation error: parent of non-root vertex %" PRId64 " is itself.\n", rank, v);
      validation_passed = 0;
    }
    if (pred[i] < -1 || pred[i] >= nglobalverts) {
      fprintf(stderr, "%d: Validation error: parent of vertex %" PRId64 " is out-of-range value %" PRId64 ".\n", rank, v, pred[i]);
      validation_passed = 0;
      any_range_errors = 1;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &any_range_errors, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

  /* Check that nvisited is correct. */
  int64_t nvisited_actual = 0;
  for (i = 0; i < nlocalverts; ++i) {
    if (pred[i] != -1) ++nvisited_actual;
  }
  MPI_Allreduce(MPI_IN_PLACE, &nvisited_actual, 1, INT64_T_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
  if (nvisited_actual != nvisited) {
    fprintf(stderr, "%d: Validation error: claimed visit count %" PRId64 " is different from actual count %" PRId64 ".\n", rank, nvisited, nvisited_actual);
    validation_passed = 0;
  }

  if (!any_range_errors) { /* Other parts of validation assume in-range values */

    /* Check that there is an edge from each vertex to its claimed
     * predecessor. */
    size_t i;
    for (i = 0; i < nlocalverts; ++i) {
      int64_t v = VERTEX_TO_GLOBAL(i);
      int64_t p = pred[i];
      if (p == -1) continue;
      int found_pred_edge = 0;
      if (v == p) found_pred_edge = 1; /* Root vertex */
      size_t ei, ei_end = g->rowstarts[i + 1];
      for (ei = g->rowstarts[i]; ei < ei_end; ++ei) {
        int64_t w = g->column[ei];
        if (w == p) {
          found_pred_edge = 1;
          break;
        }
      }
      if (!found_pred_edge) {
        fprintf(stderr, "%d: Validation error: no graph edge from vertex %" PRId64 " to its parent %" PRId64 ".\n", rank, v, p);
        validation_passed = 0;
      }
    }

    /* Create a vertex depth map to use for later validation. */
    int64_t* depth = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t));
    { /* Scope some code that has a lot of temporary variables. */
      int64_t* pred_depth = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t)); /* Depth of predecessor vertex for each local vertex */
      size_t i;
      for (i = 0; i < nlocalverts; ++i) depth[i] = INT64_MAX;
      if (root_is_mine) depth[VERTEX_LOCAL(root)] = 0;
      /* Send each vertex that appears in the local part of the predecessor map
       * to its owner; record the original locations so we can put the answers
       * into pred_depth. */
      /* Do a histogram sort by owner (this same kind of sort is used other
       * places as well).  First, count the number of vertices going to each
       * destination. */
      int* num_preds_per_owner = (int*)xcalloc(size, sizeof(int)); /* Uses zero-init */
      for (i = 0; i < nlocalverts; ++i) {
        ++num_preds_per_owner[pred[i] == -1 ? size - 1 : VERTEX_OWNER(pred[i])];
      }
      int64_t* preds_per_owner = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t)); /* Predecessors sorted by owner */
      int64_t* preds_per_owner_results_offsets = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t)); /* Indices into pred_depth to write */
      /* Second, do a prefix sum to get the displacements of the different
       * owners in the outgoing array. */
      int* pred_owner_displs = (int*)xmalloc((size + 1) * sizeof(int));
      pred_owner_displs[0] = 0;
      int r;
      for (r = 0; r < size; ++r) {
        pred_owner_displs[r + 1] = pred_owner_displs[r] + num_preds_per_owner[r];
      }
      /* Last, put the vertices into the correct positions in the array, based
       * on their owners and the counts and displacements computed earlier. */
      int* pred_owner_offsets = (int*)xmalloc((size + 1) * sizeof(int));
      memcpy(pred_owner_offsets, pred_owner_displs, (size + 1) * sizeof(int));
      for (i = 0; i < nlocalverts; ++i) {
        int* offset_ptr = &pred_owner_offsets[pred[i] == -1 ? size - 1 : VERTEX_OWNER(pred[i])];
        preds_per_owner[*offset_ptr] = pred[i];
        preds_per_owner_results_offsets[*offset_ptr] = i;
        ++*offset_ptr;
      }
      for (r = 0; r < size; ++r) {
        assert (pred_owner_offsets[r] == pred_owner_displs[r + 1]);
      }
      free(pred_owner_offsets);

      /* Send around the number of vertices that will be sent to each destination. */
      int* num_my_preds_per_sender = (int*)xmalloc(size * sizeof(int));
      MPI_Alltoall(num_preds_per_owner, 1, MPI_INT,
                   num_my_preds_per_sender, 1, MPI_INT,
                   MPI_COMM_WORLD);
      int* my_preds_per_sender_displs = (int*)xmalloc((size + 1) * sizeof(int));
      my_preds_per_sender_displs[0] = 0;
      for (r = 0; r < size; ++r) {
        my_preds_per_sender_displs[r + 1] = my_preds_per_sender_displs[r] + num_my_preds_per_sender[r];
      }
      /* Send around the actual vertex data (list of depth requests that will
       * be responded to at each BFS iteration). */
      int64_t* my_depth_requests = (int64_t*)xmalloc(my_preds_per_sender_displs[size] * sizeof(int64_t));
      int64_t* my_depth_replies = (int64_t*)xmalloc(my_preds_per_sender_displs[size] * sizeof(int64_t));
      MPI_Alltoallv(preds_per_owner, num_preds_per_owner, pred_owner_displs, INT64_T_MPI_TYPE,
                    my_depth_requests, num_my_preds_per_sender, my_preds_per_sender_displs, INT64_T_MPI_TYPE,
                    MPI_COMM_WORLD);

      int64_t* pred_depth_raw = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t)); /* Depth of predecessor vertex for each local vertex, ordered by source proc */

      /* Do a mini-BFS (naively) over just the predecessor graph (hopefully a
       * tree) produced by the real BFS; fill in the depth map. */
      while (1) {
        int any_changed = 0;
        int i;
        /* Create and send the depth values requested by other nodes.  The list
         * of requests is sent once, and are stored on the receiver so the
         * replies can be sent (possibly with updated depth values) at every
         * iteration. */
        for (i = 0; i < my_preds_per_sender_displs[size]; ++i) {
          my_depth_replies[i] = (my_depth_requests[i] == -1 ? INT64_MAX : depth[VERTEX_LOCAL(my_depth_requests[i])]);
        }
        MPI_Alltoallv(my_depth_replies, num_my_preds_per_sender, my_preds_per_sender_displs, INT64_T_MPI_TYPE,
                      pred_depth_raw, num_preds_per_owner, pred_owner_displs, INT64_T_MPI_TYPE,
                      MPI_COMM_WORLD);
        {
          size_t i;
          /* Put the received depths into the local array. */
          for (i = 0; i < nlocalverts; ++i) {
            pred_depth[preds_per_owner_results_offsets[i]] = pred_depth_raw[i];
          }
          /* Check those values to determine if they violate any correctness
           * conditions. */
          for (i = 0; i < nlocalverts; ++i) {
            int64_t v = VERTEX_TO_GLOBAL(i);
            if (v == root) {
              /* The depth and predecessor for this were checked earlier. */
            } else if (depth[i] == INT64_MAX && pred_depth[i] == INT64_MAX) {
              /* OK -- depth should be filled in later. */
            } else if (depth[i] == INT64_MAX && pred_depth[i] != INT64_MAX) {
              depth[i] = pred_depth[i] + 1;
              any_changed = 1;
            } else if (depth[i] != pred_depth[i] + 1) {
              fprintf(stderr, "%d: Validation error: BFS predecessors do not form a tree; see vertices %" PRId64 " (depth %" PRId64 ") and %" PRId64 " (depth %" PRId64 ").\n", rank, v, depth[i], pred[i], pred_depth[i]);
              validation_passed = 0;
            } else {
              /* Vertex already has its correct depth value. */
            }
          }
        }
        MPI_Allreduce(MPI_IN_PLACE, &any_changed, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        if (!any_changed) break;
      }

      free(num_preds_per_owner);
      free(num_my_preds_per_sender);
      free(preds_per_owner);
      free(preds_per_owner_results_offsets);
      free(my_preds_per_sender_displs);
      free(my_depth_requests);
      free(my_depth_replies);
      free(pred_owner_displs);
      free(pred_depth);
      free(pred_depth_raw);
    }

    /* Check that all edges connect vertices whose depths differ by at most
     * one. */
    {
      int64_t maxlocaledges = 0;
      MPI_Allreduce((void*)&nlocaledges, &maxlocaledges, 1, INT64_T_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD);
      /* We break the total list of overall edges into chunks to reduce the
       * amount of data to be sent at a time (since we are using MPI_Alltoallv
       * to send data collectively). */
      const int edge_chunk_size = (1 << 23); /* Reduce memory usage */
      int num_edge_groups = (maxlocaledges + edge_chunk_size - 1) / edge_chunk_size;
      int eg;
      for (eg = 0; eg < num_edge_groups; ++eg) {
        size_t first_edge_index = (size_t)(eg * edge_chunk_size);
        if (first_edge_index > nlocaledges) first_edge_index = nlocaledges;
        size_t last_edge_index = (size_t)((eg + 1) * edge_chunk_size);
        if (last_edge_index > nlocaledges) last_edge_index = nlocaledges;
        /* Sort the edge targets in this chunk by their owners (histogram
         * sort); see the BFS code above for details of the steps of the
         * algorithm. */
        int* num_edge_targets_by_owner = (int*)xcalloc(size, sizeof(int)); /* Uses zero-init */
        size_t ei;
        for (ei = first_edge_index; ei < last_edge_index; ++ei) {
          ++num_edge_targets_by_owner[VERTEX_OWNER(g->column[ei])];
        }
        int* edge_targets_by_owner_displs = (int*)xmalloc((size + 1) * sizeof(int));
        edge_targets_by_owner_displs[0] = 0;
        int i;
        for (i = 0; i < size; ++i) {
          edge_targets_by_owner_displs[i + 1] = edge_targets_by_owner_displs[i] + num_edge_targets_by_owner[i];
        }
        int64_t* edge_targets_by_owner = (int64_t*)xmalloc(edge_targets_by_owner_displs[size] * sizeof(int64_t));
        int64_t* edge_targets_by_owner_indices = (int64_t*)xmalloc(edge_targets_by_owner_displs[size] * sizeof(int64_t)); /* Source indices for where to write the targets */
        int* edge_targets_by_owner_offsets = (int*)xmalloc((size + 1) * sizeof(int));
        memcpy(edge_targets_by_owner_offsets, edge_targets_by_owner_displs, (size + 1) * sizeof(int));
        for (ei = first_edge_index; ei < last_edge_index; ++ei) {
          edge_targets_by_owner[edge_targets_by_owner_offsets[VERTEX_OWNER(g->column[ei])]] = g->column[ei];
          edge_targets_by_owner_indices[edge_targets_by_owner_offsets[VERTEX_OWNER(g->column[ei])]] = ei;
          ++edge_targets_by_owner_offsets[VERTEX_OWNER(g->column[ei])];
        }
        for (i = 0; i < size; ++i) {
          assert (edge_targets_by_owner_offsets[i] == edge_targets_by_owner_displs[i + 1]);
        }
        free(edge_targets_by_owner_offsets);

        /* Send around the number of data elements that will be sent later. */
        int* num_incoming_targets_by_src = (int*)xmalloc(size * sizeof(int));
        MPI_Alltoall(num_edge_targets_by_owner, 1, MPI_INT,
                     num_incoming_targets_by_src, 1, MPI_INT,
                     MPI_COMM_WORLD);
        int* incoming_targets_by_src_displs = (int*)xmalloc((size + 1) * sizeof(int));
        incoming_targets_by_src_displs[0] = 0;
        for (i = 0; i < size; ++i) {
          incoming_targets_by_src_displs[i + 1] = incoming_targets_by_src_displs[i] + num_incoming_targets_by_src[i];
        }

        int64_t* target_depth_requests = (int64_t*)xmalloc(incoming_targets_by_src_displs[size] * sizeof(int64_t));
        int64_t* target_depth_replies = (int64_t*)xmalloc(incoming_targets_by_src_displs[size] * sizeof(int64_t));

        /* Send the actual requests for the depths of edge targets. */
        MPI_Alltoallv(edge_targets_by_owner, num_edge_targets_by_owner, edge_targets_by_owner_displs, INT64_T_MPI_TYPE,
                      target_depth_requests, num_incoming_targets_by_src, incoming_targets_by_src_displs, INT64_T_MPI_TYPE,
                      MPI_COMM_WORLD);

        free(edge_targets_by_owner);

        /* Fill in the replies for the requests sent to me. */
        for (i = 0; i < incoming_targets_by_src_displs[size]; ++i) {
          assert (VERTEX_OWNER(target_depth_requests[i]) == rank);
          target_depth_replies[i] = depth[VERTEX_LOCAL(target_depth_requests[i])];
        }

        free(target_depth_requests);

        int64_t* target_depth_raw = (int64_t*)xmalloc((last_edge_index - first_edge_index) * sizeof(int64_t));

        /* Send back the replies. */
        MPI_Alltoallv(target_depth_replies, num_incoming_targets_by_src, incoming_targets_by_src_displs, INT64_T_MPI_TYPE,
                      target_depth_raw, num_edge_targets_by_owner, edge_targets_by_owner_displs, INT64_T_MPI_TYPE,
                      MPI_COMM_WORLD);

        free(target_depth_replies);
        free(num_incoming_targets_by_src);
        free(num_edge_targets_by_owner);
        free(incoming_targets_by_src_displs);
        free(edge_targets_by_owner_displs);

        int64_t* target_depth = (int64_t*)xmalloc((last_edge_index - first_edge_index) * sizeof(int64_t));

        /* Put the replies into the proper order (original order of the edges).
         * */
        for (ei = 0; ei < last_edge_index - first_edge_index; ++ei) {
          target_depth[edge_targets_by_owner_indices[ei] - first_edge_index] = target_depth_raw[ei];
        }

        free(target_depth_raw);
        free(edge_targets_by_owner_indices);

        /* Check the depth relationship of the endpoints of each edge in the
         * current chunk. */
        size_t src_i = 0;
        for (ei = first_edge_index; ei < last_edge_index; ++ei) {
          while (ei >= g->rowstarts[src_i + 1]) {
            ++src_i;
          }
          int64_t src = VERTEX_TO_GLOBAL(src_i);
          int64_t src_depth = depth[src_i];
          int64_t tgt = g->column[ei];
          int64_t tgt_depth = target_depth[ei - first_edge_index];
          if (src_depth != INT64_MAX && tgt_depth == INT64_MAX) {
            fprintf(stderr, "%d: Validation error: edge connects vertex %" PRId64 " in the BFS tree (depth %" PRId64 ") to vertex %" PRId64 " outside the tree.\n", rank, src, src_depth, tgt);
            validation_passed = 0;
          } else if (src_depth == INT64_MAX && tgt_depth != INT64_MAX) {
            /* Skip this for now; this problem will be caught when scanning
             * reversed copy of this edge.  Set the failure flag, though,
             * just in case. */
            validation_passed = 0;
          } else if (src_depth - tgt_depth < -1 ||
                     src_depth - tgt_depth > 1) {
            fprintf(stderr, "%d: Validation error: depths of edge endpoints %" PRId64 " (depth %" PRId64 ") and %" PRId64 " (depth %" PRId64 ") are too far apart (abs. val. > 1).\n", rank, src, src_depth, tgt, tgt_depth);
            validation_passed = 0;
          }
        }
        free(target_depth);
      }
    }

    free(depth);

  } /* End of part skipped by range errors */
  
  /* Collect the global validation result. */
  MPI_Allreduce(MPI_IN_PLACE, &validation_passed, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  return validation_passed;
}
