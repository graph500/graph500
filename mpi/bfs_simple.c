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

/* This version is the traditional level-synchronized BFS using two queues.  A
 * bitmap is used to indicate which vertices have been visited.  Messages are
 * sent and processed asynchronously throughout the code to hopefully overlap
 * communication with computation. */
void run_mpi_bfs(const csr_graph* const g, int64_t root, int64_t* pred, int64_t* nvisited) {
  const size_t nlocalverts = g->nlocalverts;
  int64_t nvisited_local = 0;

  /* Set up the queues. */
  int64_t* oldq = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t));
  int64_t* newq = (int64_t*)xmalloc(nlocalverts * sizeof(int64_t));
  size_t oldq_count = 0;
  size_t newq_count = 0;

  /* Set up the visited bitmap. */
  const int ulong_bits = sizeof(unsigned long) * CHAR_BIT;
  int64_t visited_size = (nlocalverts + ulong_bits - 1) / ulong_bits;
  unsigned long* visited = (unsigned long*)xcalloc(visited_size, sizeof(unsigned long)); /* Uses zero-init */
#define SET_VISITED(v) do {visited[VERTEX_LOCAL((v)) / ulong_bits] |= (1UL << (VERTEX_LOCAL((v)) % ulong_bits));} while (0)
#define TEST_VISITED(v) ((visited[VERTEX_LOCAL((v)) / ulong_bits] & (1UL << (VERTEX_LOCAL((v)) % ulong_bits))) != 0)

  /* Set up buffers for message coalescing, MPI requests, etc. for
   * communication. */
  const int coalescing_size = 256;
  int64_t* outgoing = (int64_t*)xMPI_Alloc_mem(coalescing_size * size * 2 * sizeof(int64_t));
  size_t* outgoing_counts = (size_t*)xmalloc(size * sizeof(size_t)) /* 2x actual count */;
  MPI_Request* outgoing_reqs = (MPI_Request*)xmalloc(size * sizeof(MPI_Request));
  int* outgoing_reqs_active = (int*)xcalloc(size, sizeof(int)); /* Uses zero-init */
  int64_t* recvbuf = (int64_t*)xMPI_Alloc_mem(coalescing_size * 2 * sizeof(int64_t));
  MPI_Request recvreq;
  int recvreq_active = 0;

  /* Termination counter for each level: this variable counts the number of
   * ranks that have said that they are done sending to me in the current
   * level.  This rank can stop listening for new messages when it reaches
   * size. */
  int num_ranks_done;

  /* Set all vertices to "not visited." */
  {size_t i; for (i = 0; i < nlocalverts; ++i) pred[i] = -1;}

  /* Mark the root and put it into the queue. */
  if (VERTEX_OWNER(root) == rank) {
    SET_VISITED(root);
    pred[VERTEX_LOCAL(root)] = root;
    oldq[oldq_count++] = root;
  }

#define CHECK_MPI_REQS \
  /* Check all MPI requests and handle any that have completed. */ \
  do { \
    /* Test for incoming vertices to put onto the queue. */ \
    while (recvreq_active) { \
      int flag; \
      MPI_Status st; \
      MPI_Test(&recvreq, &flag, &st); \
      if (flag) { \
        recvreq_active = 0; \
        int count; \
        MPI_Get_count(&st, INT64_T_MPI_TYPE, &count); \
        /* count == 0 is a signal from a rank that it is done sending to me
         * (using MPI's non-overtaking rules to keep that signal after all
         * "real" messages. */ \
        if (count == 0) { \
          ++num_ranks_done; \
        } else { \
          int j; \
          for (j = 0; j < count; j += 2) { \
            int64_t tgt = recvbuf[j]; \
            int64_t src = recvbuf[j + 1]; \
            /* Process one incoming edge. */ \
            assert (VERTEX_OWNER(tgt) == rank); \
            if (!TEST_VISITED(tgt)) { \
              SET_VISITED(tgt); \
              pred[VERTEX_LOCAL(tgt)] = src; \
              newq[newq_count++] = tgt; \
            } \
          } \
        } \
        /* Restart the receive if more messages will be coming. */ \
        if (num_ranks_done < size) { \
          MPI_Irecv(recvbuf, coalescing_size * 2, INT64_T_MPI_TYPE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recvreq); \
          recvreq_active = 1; \
        } \
      } else break; \
    } \
    /* Mark any sends that completed as inactive so their buffers can be
     * reused. */ \
    int c; \
    for (c = 0; c < size; ++c) { \
      if (outgoing_reqs_active[c]) { \
        int flag; \
        MPI_Test(&outgoing_reqs[c], &flag, MPI_STATUS_IGNORE); \
        if (flag) outgoing_reqs_active[c] = 0; \
      } \
    } \
  } while (0)

  while (1) {
    memset(outgoing_counts, 0, size * sizeof(size_t));
    num_ranks_done = 1; /* I never send to myself, so I'm always done */
    
    /* Start the initial receive. */
    if (num_ranks_done < size) {
      MPI_Irecv(recvbuf, coalescing_size * 2, INT64_T_MPI_TYPE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recvreq);
      recvreq_active = 1;
    }

    /* Step through the current level's queue. */
    size_t i;
    for (i = 0; i < oldq_count; ++i) {
      CHECK_MPI_REQS;
      assert (VERTEX_OWNER(oldq[i]) == rank);
      assert (pred[VERTEX_LOCAL(oldq[i])] >= 0 && pred[VERTEX_LOCAL(oldq[i])] < g->nglobalverts);
      ++nvisited_local;
      int64_t src = oldq[i];
      /* Iterate through its incident edges. */
      size_t j, j_end = g->rowstarts[VERTEX_LOCAL(oldq[i]) + 1];
      for (j = g->rowstarts[VERTEX_LOCAL(oldq[i])]; j < j_end; ++j) {
        int64_t tgt = g->column[j];
        int owner = VERTEX_OWNER(tgt);
        /* If the other endpoint is mine, update the visited map, predecessor
         * map, and next-level queue locally; otherwise, send the target and
         * the current vertex (its possible predecessor) to the target's owner.
         * */
        if (owner == rank) {
          if (!TEST_VISITED(tgt)) {
            SET_VISITED(tgt);
            pred[VERTEX_LOCAL(tgt)] = src;
            newq[newq_count++] = tgt;
          }
        } else {
          while (outgoing_reqs_active[owner]) CHECK_MPI_REQS; /* Wait for buffer to be available */
          size_t c = outgoing_counts[owner];
          outgoing[owner * coalescing_size * 2 + c] = tgt;
          outgoing[owner * coalescing_size * 2 + c + 1] = src;
          outgoing_counts[owner] += 2;
          if (outgoing_counts[owner] == coalescing_size * 2) {
            MPI_Isend(&outgoing[owner * coalescing_size * 2], coalescing_size * 2, INT64_T_MPI_TYPE, owner, 0, MPI_COMM_WORLD, &outgoing_reqs[owner]);
            outgoing_reqs_active[owner] = 1;
            outgoing_counts[owner] = 0;
          }
        }
      }
    }
    /* Flush any coalescing buffers that still have messages. */
    int offset;
    for (offset = 1; offset < size; ++offset) {
      int dest = MOD_SIZE(rank + offset);
      if (outgoing_counts[dest] != 0) {
        while (outgoing_reqs_active[dest]) CHECK_MPI_REQS;
        MPI_Isend(&outgoing[dest * coalescing_size * 2], outgoing_counts[dest], INT64_T_MPI_TYPE, dest, 0, MPI_COMM_WORLD, &outgoing_reqs[dest]);
        outgoing_reqs_active[dest] = 1;
        outgoing_counts[dest] = 0;
      }
      /* Wait until all sends to this destination are done. */
      while (outgoing_reqs_active[dest]) CHECK_MPI_REQS;
      /* Tell the destination that we are done sending to them. */
      MPI_Isend(&outgoing[dest * coalescing_size * 2], 0, INT64_T_MPI_TYPE, dest, 0, MPI_COMM_WORLD, &outgoing_reqs[dest]); /* Signal no more sends */
      outgoing_reqs_active[dest] = 1;
      while (outgoing_reqs_active[dest]) CHECK_MPI_REQS;
    }
    /* Wait until everyone else is done (and thus couldn't send us any more
     * messages). */
    while (num_ranks_done < size) CHECK_MPI_REQS;

    /* Test globally if all queues are empty. */
    int64_t global_newq_count;
    MPI_Allreduce(&newq_count, &global_newq_count, 1, INT64_T_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);

    /* Quit if they all are empty. */
    if (global_newq_count == 0) break;

    /* Swap old and new queues; clear new queue for next level. */
    {int64_t* temp = oldq; oldq = newq; newq = temp;}
    oldq_count = newq_count;
    newq_count = 0;
  }
#undef CHECK_MPI_REQS

  /* Add up total visited-vertex count. */
  MPI_Allreduce(MPI_IN_PLACE, &nvisited_local, 1, INT64_T_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
  *nvisited = nvisited_local;

  free(oldq);
  free(newq);
  MPI_Free_mem(outgoing);
  free(outgoing_counts);
  free(outgoing_reqs);
  free(outgoing_reqs_active);
  MPI_Free_mem(recvbuf);
  free(visited);
}
