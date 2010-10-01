/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

/* MPI code to apply a distributed permutation.  The versions for sequential,
 * XMT, and OpenMP are trivial and so are in the corresponding versions of
 * make_graph() (make_graph.c), but the MPI one is long and so it is in a
 * separate file. */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#include "graph_generator.h"
#include "permutation_gen.h"
#include "apply_permutation_mpi.h"
#include "utils.h"

void gather_block_distribution_info(MPI_Comm comm, const int64_t local_perm_size, const int64_t global_perm_size, int64_t* perm_displs /* size = MPI comm size + 1 */, int** perm_owner_table_ptr /* malloc'ed in here */, int64_t** perm_owner_cutoff_ptr /* malloc'ed in here */, int* lg_minpermsize_ptr, int64_t* maxpermsize_ptr) {
  int rank, size;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  /* Get the size of the permutation fragment on each node. */
  int64_t* perm_sizes = (int64_t*)xmalloc(size * sizeof(int64_t));
  MPI_Allgather((void*)&local_perm_size, 1, INT64_T_MPI_TYPE, perm_sizes, 1, INT64_T_MPI_TYPE, comm);

  /* Sanity check the permutation to make sure its distribution is reasonable
   * (PRNG problems can cause bad distributions) and get some statistics to use
   * later. */
  int minpermsize = local_perm_size, maxpermsize = local_perm_size;
  int64_t i;
  for (i = 0; i < size; ++i) {
    if (perm_sizes[i] < minpermsize) minpermsize = perm_sizes[i];
    if (perm_sizes[i] > maxpermsize) maxpermsize = perm_sizes[i];
  }
  *maxpermsize_ptr = maxpermsize;

  if (minpermsize == 0) {
    fprintf(stderr, "One node has permutation size 0\n");
    abort();
  }
  if (maxpermsize > 2 * minpermsize && rank == 0 && maxpermsize >= 1000) {
    fprintf(stderr, "Minimum permutation bucket size %d is too small compared to maximum size %d\n", minpermsize, maxpermsize);
    abort();
  }

  /* Compute floor(log_2(minpermsize)) for fast computation in lookup table. */
  int lg_minpermsize = 0;
  while ((1 << lg_minpermsize) <= minpermsize) ++lg_minpermsize;
  --lg_minpermsize;
  int minpermsize_for_table = (1 << lg_minpermsize);
  *lg_minpermsize_ptr = lg_minpermsize;

  /* Prefix sum permutation sizes (exclusive) for owner lookups. */
  perm_displs[0] = 0;
  for (i = 1; i < size + 1; ++i) perm_displs[i] = perm_displs[i - 1] + perm_sizes[i - 1];
  free(perm_sizes); perm_sizes = NULL;

  /* Create tables containing permutation_owner(minpermsize_for_table * i) and
   * perm_displs[permutation_owner(minpermsize_for_table * i) + 1] for i = 0 ...  global_perm_size /
   * minpermsize_for_table.  These two tables allow a constant-time owner lookup. */
  size_t perm_owner_table_size = (size_t)((global_perm_size + minpermsize_for_table - 1) >> lg_minpermsize);
  int* perm_owner_table = (int*)xmalloc(perm_owner_table_size * sizeof(int));
  *perm_owner_table_ptr = perm_owner_table;
  int64_t* perm_owner_cutoff = (int64_t*)xmalloc(perm_owner_table_size * sizeof(int64_t));
  *perm_owner_cutoff_ptr = perm_owner_cutoff;
  int cur_owner = 0;
  for (i = 0; i < perm_owner_table_size; ++i) {
    assert (cur_owner < size);
    while ((i << lg_minpermsize) >= perm_displs[cur_owner + 1]) {
      ++cur_owner;
      assert (cur_owner < size);
    }
    perm_owner_table[i] = cur_owner;
    perm_owner_cutoff[i] = perm_displs[cur_owner + 1];
    assert (perm_displs[cur_owner] <= (i << lg_minpermsize) &&
            perm_displs[cur_owner + 1] > (i << lg_minpermsize));
  }

}

void apply_permutation_mpi(MPI_Comm comm, const int64_t local_perm_size, const int64_t* const local_vertex_perm, const int64_t N, const int64_t nedges, int64_t* result) {
  int rank, size;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int64_t* perm_displs = (int64_t*)xmalloc((size + 1) * sizeof(int64_t));
  int* perm_owner_table;
  int64_t* perm_owner_cutoff;
  int lg_minpermsize;
  int64_t maxpermsize;
  gather_block_distribution_info(comm, local_perm_size, N, perm_displs, &perm_owner_table, &perm_owner_cutoff, &lg_minpermsize, &maxpermsize);

  int64_t my_perm_displ = perm_displs[rank];

#define LOOKUP_PERM_OWNER(v) \
  (perm_owner_table[(v) >> lg_minpermsize] + \
   ((v) >= perm_owner_cutoff[(v) >> lg_minpermsize]))

  /* From here on until a cleanup at the end of permutation application, vertex
   * numbers in edges will be coded as either (unpermuted vertex) or (permuted
   * vertex - N - 1) [avoiding -1 in the second case because that still means
   * unused vertex slot]. */

  /* Apply any permutation elements that I own. */
  int64_t i;
  for (i = 0; i < 2 * nedges; ++i) {
    if (result[i] >= perm_displs[rank] && result[i] < perm_displs[rank + 1]) {
      result[i] = local_vertex_perm[result[i] - perm_displs[rank]] - N - 1;
    }
  }

  /* Create a bitmap (as bytes) for which nodes I will need data from. */
  unsigned char* owner_bitmap = (unsigned char*)xcalloc(size, sizeof(unsigned char)); /* Uses zero-init */
  for (i = 0; i < 2 * nedges; ++i) {
    if (result[i] >= 0) {
      owner_bitmap[LOOKUP_PERM_OWNER(result[i])] = 1;
    }
  }

  /* Send those bitmaps so owners of permutation segments know who they will be
   * sending data to. */
  unsigned char* send_bitmap = (unsigned char*)xmalloc(size * sizeof(unsigned char));
  MPI_Alltoall(owner_bitmap, 1, MPI_UNSIGNED_CHAR,
               send_bitmap, 1, MPI_UNSIGNED_CHAR,
               comm);

  /* Go through node offsets in batches of node_chunk_size.  Create bitmaps of
   * which particular permutation elements will be needed from each of those
   * nodes, and which local vertices in result should be checked when those
   * elements arrive. */
  const int node_chunk_size = 64;
  const int ulong_bits = sizeof(unsigned long) * CHAR_BIT;
  const size_t needed_element_bitmap_size = (maxpermsize + ulong_bits - 1) / ulong_bits; /* Per destination */
  unsigned long* needed_element_bitmap =
    (unsigned long*)xmalloc(needed_element_bitmap_size * node_chunk_size * sizeof(unsigned long));
  /* Number of needed elements per source */
  int64_t* num_needed_elements =
    (int64_t*)xmalloc(node_chunk_size * sizeof(int64_t));
  const size_t rescan_bitmap_size = (2 * nedges + ulong_bits - 1) / ulong_bits;
  unsigned long* rescan_bitmap =
    (unsigned long*)xmalloc(rescan_bitmap_size * sizeof(unsigned long));
  int64_t* incoming_perm_requests =
    (int64_t*)xmalloc(maxpermsize * sizeof(int64_t));
  int64_t* incoming_perm_replies =
    (int64_t*)xmalloc(maxpermsize * sizeof(int64_t));
  int64_t* outgoing_perm_requests =
    (int64_t*)xmalloc(maxpermsize * sizeof(int64_t));
  int64_t* outgoing_perm_replies =
    (int64_t*)xmalloc(maxpermsize * sizeof(int64_t));
  int64_t* perm_replies_spread =
    (int64_t*)xmalloc(maxpermsize * sizeof(int64_t)); /* Replies scattered for direct indexing */

  int offset;
  for (offset = 1 /* Already did our own perm elements */;
       offset < size; offset += node_chunk_size) {
    /* Chunk of destinations with this offset is [rank+offset, rank+min(size,
     * offset+node_chunk_size)), all of those modulo size. */
    int num_nodes_being_processed = size - offset;
    if (num_nodes_being_processed > node_chunk_size) {
      num_nodes_being_processed = node_chunk_size;
    }

    /* Clear bitmaps. */
    memset(num_needed_elements, 0, node_chunk_size * sizeof(int64_t));
    memset(needed_element_bitmap, 0, node_chunk_size * needed_element_bitmap_size * sizeof(unsigned long));
    memset(rescan_bitmap, 0, rescan_bitmap_size * sizeof(unsigned long));
#define SET_BIT(bm, idx) do { \
                           (bm)[(idx) / ulong_bits] |= (1UL << ((idx) % ulong_bits)); \
                         } while (0)
#define TEST_BIT_SCALAR(bm, idx) (((bm) & (1UL << ((idx) % ulong_bits))) != 0)
#define TEST_BIT(bm, idx) TEST_BIT_SCALAR((bm)[(idx) / ulong_bits], idx)

    /* Fill in needed_element_bitmap (and num_needed_elements), as well as
     * rescan_bitmap, for this set of destinations (permutation owners). */
    for (i = 0; i < 2 * nedges; ++i) {
      int64_t v = result[i];
      if (v < 0) continue;
      int owner = LOOKUP_PERM_OWNER(v);
      int64_t owner_displ = perm_displs[owner];
      int delta = (owner + size - rank) % size; /* Rank difference from mine, modulo size */
      if (delta >= offset && delta < offset + node_chunk_size) {
        unsigned long* bm = needed_element_bitmap + (delta - offset) * needed_element_bitmap_size;
        if (!TEST_BIT(bm, v - owner_displ)) {
          SET_BIT(bm, v - owner_displ);
          ++num_needed_elements[delta - offset];
        }
        SET_BIT(rescan_bitmap, i);
      }
    }

    /* For each delta in the current range, set up the sends and receives if
     * necessary (as determined by send_bitmap and owner_bitmap). */
    int delta;
    for (delta = offset; delta < offset + num_nodes_being_processed; ++delta) {
      MPI_Request incomingreq;
      int src = (rank + size - delta) % size;
      if (send_bitmap[src]) { /* Will I get any requests from src? */
        /* Set up the receive for incoming permutation data requests. */
        MPI_Irecv(incoming_perm_requests, maxpermsize, INT64_T_MPI_TYPE, src, 0, comm, &incomingreq);
      }
      int dest = (rank + delta) % size;
      MPI_Request responsereq;
      int output_offset = 0;
      int64_t perm_displ_for_dest = perm_displs[dest]; /* First permutation element on dest */
      if (owner_bitmap[dest]) {
        /* Convert the needed_element_bitmap chunk for this destination into an
         * array of individual vertex indices. */
        const unsigned long* bm = needed_element_bitmap + (delta - offset) * needed_element_bitmap_size;
        int64_t block;
        for (block = 0; block < needed_element_bitmap_size; ++block) {
          unsigned long val = bm[block];
          if (!val) continue;
          int bit;
          for (bit = 0; bit < ulong_bits; ++bit) {
            if (!TEST_BIT_SCALAR(val, bit)) continue;
            outgoing_perm_requests[output_offset++] = perm_displ_for_dest + block * ulong_bits + bit;
          }
        }
        /* Do irecv first to allow rsend on destination. */
        MPI_Irecv(outgoing_perm_replies, output_offset, INT64_T_MPI_TYPE, dest, 1, comm, &responsereq);
        /* Send the requests, using ssend for flow control (this can be changed
         * if your system has good flow control). */
        MPI_Ssend(outgoing_perm_requests, output_offset, INT64_T_MPI_TYPE, dest, 0, comm);
      }
      if (send_bitmap[src]) {
        /* Wait for and process an incoming request if send_bitmap says there
         * will be one from this src. */
        MPI_Status st;
        MPI_Wait(&incomingreq, &st);
        int count;
        MPI_Get_count(&st, INT64_T_MPI_TYPE, &count);
        int i;
        for (i = 0; i < count; ++i) {
          incoming_perm_replies[i] = local_vertex_perm[incoming_perm_requests[i] - my_perm_displ];
        }
        /* Rsend the reply for performance; note recv being before send in
         * request-sending code above. */
        MPI_Rsend(incoming_perm_replies, count, INT64_T_MPI_TYPE, src, 1, comm);
      }
      if (owner_bitmap[dest]) {
        MPI_Wait(&responsereq, MPI_STATUS_IGNORE);
        /* Distribute the (request, reply) pairs into a large array for easy
         * random access without searching. */
        int i;
        for (i = 0; i < output_offset; ++i) {
          perm_replies_spread[outgoing_perm_requests[i] - perm_displ_for_dest] = outgoing_perm_replies[i];
        }
        /* For those edges that are marked as needing rescanning for this group
         * of destinations, update any with endpoints that belong to the
         * current destination. */
        int64_t block;
        for (block = 0; block < rescan_bitmap_size; ++block) {
          unsigned long val = rescan_bitmap[block];
          if (val == 0) continue;
          int bit;
          for (bit = 0; bit < ulong_bits; ++bit) {
            if (!TEST_BIT_SCALAR(val, bit)) continue;
            int64_t i = block * ulong_bits + bit;
            int64_t v = result[i];
            if (v >= 0 && LOOKUP_PERM_OWNER(v) == dest) {
              result[i] = perm_replies_spread[v - perm_displ_for_dest] - N - 1;
            }
          }
        }
      }
    }
#undef SET_BIT
#undef TEST_BIT_SCALAR
#undef TEST_BIT
#undef LOOKUP_PERM_OWNER
  }
  free(perm_owner_table); perm_owner_table = NULL;
  free(perm_owner_cutoff); perm_owner_cutoff = NULL;
  free(perm_displs); perm_displs = NULL;
  free(rescan_bitmap); rescan_bitmap = NULL;
  free(needed_element_bitmap); needed_element_bitmap = NULL;
  free(num_needed_elements); num_needed_elements = NULL;
  free(incoming_perm_requests); incoming_perm_requests = NULL;
  free(incoming_perm_replies); incoming_perm_replies = NULL;
  free(outgoing_perm_requests); outgoing_perm_requests = NULL;
  free(outgoing_perm_replies); outgoing_perm_replies = NULL;
  free(perm_replies_spread); perm_replies_spread = NULL;

  /* Undo bias of permuted vertices in result */
  for (i = 0; i < 2 * nedges; i += 2) {
    if (result[i] != -1) {
      int64_t v1 = result[i] + N + 1;
      int64_t v2 = result[i + 1] + N + 1;
      assert (v1 >= 0 && v1 < N);
      assert (v2 >= 0 && v2 < N);
      /* Sort these since otherwise the directions of the permuted edges would
       * give away the unscrambled vertex order. */
      result[i] = (v1 < v2) ? v1 : v2;
      result[i + 1] = (v1 < v2) ? v2 : v1;
    }
  }

  free(send_bitmap); send_bitmap = NULL;
  free(owner_bitmap); owner_bitmap = NULL;
}
