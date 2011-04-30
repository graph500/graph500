/* Copyright (C) 2010 The Trustees of Indiana University.                  */
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

static oned_csr_graph g;

void make_graph_data_structure(const tuple_graph* const tg) {
  convert_graph_to_oned_csr(tg, &g);
}

void free_graph_data_structure(void) {
  free_oned_csr_graph(&g);
}

int bfs_writes_depth_map(void) {
  return 0;
}

/* This BFS represents its queues as bitmaps and uses some data representation
 * tricks to fit with the use of MPI one-sided operations.  It is not much
 * faster than the standard version on the machines I have tested it on, but
 * systems that have good RDMA hardware and good MPI one-sided implementations
 * might get better performance from it.  This code might also be good to
 * translate to UPC, Co-array Fortran, SHMEM, or GASNet since those systems are
 * more designed for one-sided remote memory operations. */
void run_bfs(int64_t root, int64_t* pred) {
  const size_t nlocalverts = g.nlocalverts;
  const int64_t nglobalverts = g.nglobalverts;

  /* Set up a second predecessor map so we can read from one and modify the
   * other. */
  int64_t* orig_pred = pred;
  int64_t* pred2 = (int64_t*)xMPI_Alloc_mem(nlocalverts * sizeof(int64_t));

  /* The queues (old and new) are represented as bitmaps.  Each bit in the
   * queue bitmap says to check elts_per_queue_bit elements in the predecessor
   * map for vertices that need to be visited.  In other words, the queue
   * bitmap is an overapproximation of the actual queue; because MPI_Accumulate
   * does not get any information on the result of the update, sometimes
   * elements are also added to the bitmap when they were actually already
   * black.  Because of this, the predecessor map needs to be checked to be
   * sure a given vertex actually needs to be processed. */
  const int elts_per_queue_bit = 4;
  const int ulong_bits = sizeof(unsigned long) * CHAR_BIT;
  int64_t queue_nbits = (nlocalverts + elts_per_queue_bit - 1) / elts_per_queue_bit;
  int64_t queue_nwords = (queue_nbits + ulong_bits - 1) / ulong_bits;
  unsigned long* queue_bitmap1 = (unsigned long*)xMPI_Alloc_mem(queue_nwords * sizeof(unsigned long));
  unsigned long* queue_bitmap2 = (unsigned long*)xMPI_Alloc_mem(queue_nwords * sizeof(unsigned long));
  memset(queue_bitmap1, 0, queue_nwords * sizeof(unsigned long));

  /* List of local vertices (used as sources in MPI_Accumulate). */
  int64_t* local_vertices = (int64_t*)xMPI_Alloc_mem(nlocalverts * sizeof(int64_t));
  {size_t i; for (i = 0; i < nlocalverts; ++i) local_vertices[i] = VERTEX_TO_GLOBAL(rank, i);}

  /* List of all bit masks for an unsigned long (used as sources in
   * MPI_Accumulate). */
  unsigned long masks[ulong_bits];
  {int i; for (i = 0; i < ulong_bits; ++i) masks[i] = (1UL << i);}

  /* Coding of predecessor map: */
  /* - White (not visited): INT64_MAX */
  /* - Grey (in queue): 0 .. nglobalverts-1 */
  /* - Black (done): -nglobalverts .. -1 */

  /* Set initial predecessor map. */
  {size_t i; for (i = 0; i < nlocalverts; ++i) pred[i] = INT64_MAX;}

  /* Mark root as grey and add it to the queue. */
  if (VERTEX_OWNER(root) == rank) {
    pred[VERTEX_LOCAL(root)] = root;
    queue_bitmap1[VERTEX_LOCAL(root) / elts_per_queue_bit / ulong_bits] |= (1UL << ((VERTEX_LOCAL(root) / elts_per_queue_bit) % ulong_bits));
  }

  /* Create MPI windows on the two predecessor arrays and the two queues. */
  MPI_Win pred_win, pred2_win, queue1_win, queue2_win;
  MPI_Win_create(pred, nlocalverts * sizeof(int64_t), sizeof(int64_t), MPI_INFO_NULL, MPI_COMM_WORLD, &pred_win);
  MPI_Win_create(pred2, nlocalverts * sizeof(int64_t), sizeof(int64_t), MPI_INFO_NULL, MPI_COMM_WORLD, &pred2_win);
  MPI_Win_create(queue_bitmap1, queue_nwords * sizeof(unsigned long), sizeof(unsigned long), MPI_INFO_NULL, MPI_COMM_WORLD, &queue1_win);
  MPI_Win_create(queue_bitmap2, queue_nwords * sizeof(unsigned long), sizeof(unsigned long), MPI_INFO_NULL, MPI_COMM_WORLD, &queue2_win);

  while (1) {
    int64_t i;
    /* Clear the next-level queue. */
    memset(queue_bitmap2, 0, queue_nwords * sizeof(unsigned long));

    /* The pred2 array is pred with all grey vertices changed to black. */
    memcpy(pred2, pred, nlocalverts * sizeof(int64_t));
    for (i = 0; i < (int64_t)nlocalverts; ++i) {
      if (pred2[i] >= 0 && pred2[i] < nglobalverts) pred2[i] -= nglobalverts;
    }

    /* Start one-sided operations for this level. */
    MPI_Win_fence(MPI_MODE_NOPRECEDE, pred2_win);
    MPI_Win_fence(MPI_MODE_NOPRECEDE, queue2_win);

    /* Step through the words of the queue bitmap. */
    for (i = 0; i < queue_nwords; ++i) {
      unsigned long val = queue_bitmap1[i];
      int bitnum;
      /* Skip any that are all zero. */
      if (!val) continue;
      /* Scan the bits in the word. */
      for (bitnum = 0; bitnum < ulong_bits; ++bitnum) {
        size_t first_v_local = (size_t)((i * ulong_bits + bitnum) * elts_per_queue_bit);
        if (first_v_local >= nlocalverts) break;
        int bit = (int)((val >> bitnum) & 1);
        /* Skip any that are zero. */
        if (!bit) continue;
        /* Scan the queue elements corresponding to this bit. */
        int qelem_idx;
        for (qelem_idx = 0; qelem_idx < elts_per_queue_bit; ++qelem_idx) {
          size_t v_local = first_v_local + qelem_idx;
          if (v_local >= nlocalverts) continue;
          /* Since the queue is an overapproximation, check the predecessor map
           * to be sure this vertex is grey. */
          if (pred[v_local] >= 0 && pred[v_local] < nglobalverts) {
            size_t ei, ei_end = g.rowstarts[v_local + 1];
            /* Walk the incident edges. */
            for (ei = g.rowstarts[v_local]; ei < ei_end; ++ei) {
              int64_t w = g.column[ei];
              if (w == VERTEX_TO_GLOBAL(rank, v_local)) continue; /* Self-loop */
              /* Set the predecessor of the other edge endpoint (note use of
               * MPI_MIN and the coding of the predecessor map). */
              MPI_Accumulate(&local_vertices[v_local], 1, MPI_INT64_T, VERTEX_OWNER(w), VERTEX_LOCAL(w), 1, MPI_INT64_T, MPI_MIN, pred2_win);
              /* Mark the endpoint in the remote queue (note that the min may
               * not do an update, so the queue is an overapproximation in this
               * way as well). */
              MPI_Accumulate(&masks[((VERTEX_LOCAL(w) / elts_per_queue_bit) % ulong_bits)], 1, MPI_UNSIGNED_LONG, VERTEX_OWNER(w), VERTEX_LOCAL(w) / elts_per_queue_bit / ulong_bits, 1, MPI_UNSIGNED_LONG, MPI_BOR, queue2_win);
            }
          }
        }
      }
    }
    /* End one-sided operations. */
    MPI_Win_fence(MPI_MODE_NOSUCCEED, queue2_win);
    MPI_Win_fence(MPI_MODE_NOSUCCEED, pred2_win);

    /* Test if there are any elements in the next-level queue (globally); stop
     * if none. */
    int any_set = 0;
    for (i = 0; i < queue_nwords; ++i) {
      if (queue_bitmap2[i] != 0) {any_set = 1; break;}
    }
    MPI_Allreduce(MPI_IN_PLACE, &any_set, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    if (!any_set) break;

    /* Swap queues and predecessor maps. */
    {MPI_Win temp = queue1_win; queue1_win = queue2_win; queue2_win = temp;}
    {unsigned long* temp = queue_bitmap1; queue_bitmap1 = queue_bitmap2; queue_bitmap2 = temp;}
    {MPI_Win temp = pred_win; pred_win = pred2_win; pred2_win = temp;}
    {int64_t* temp = pred; pred = pred2; pred2 = temp;}
  }
  MPI_Win_free(&pred_win);
  MPI_Win_free(&pred2_win);
  MPI_Win_free(&queue1_win);
  MPI_Win_free(&queue2_win);
  MPI_Free_mem(local_vertices);
  MPI_Free_mem(queue_bitmap1);
  MPI_Free_mem(queue_bitmap2);

  /* Clean up the predecessor map swapping since the surrounding code does not
   * allow the BFS to change the predecessor map pointer. */
  if (pred2 != orig_pred) {
    memcpy(orig_pred, pred2, nlocalverts * sizeof(int64_t));
    MPI_Free_mem(pred2);
  } else {
    MPI_Free_mem(pred);
  }

  /* Change from special coding of predecessor map to the one the benchmark
   * requires. */
  size_t i;
  for (i = 0; i < nlocalverts; ++i) {
    if (orig_pred[i] < 0) {
      orig_pred[i] += nglobalverts;
    } else if (orig_pred[i] == INT64_MAX) {
      orig_pred[i] = -1;
    }
  }
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
