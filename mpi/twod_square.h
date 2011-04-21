/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef TWOD_SQUARE_H
#define TWOD_SQUARE_H

#include "common.h"
#include <limits.h>
#include <assert.h>

#define ULONG_BITS (sizeof(unsigned long) * CHAR_BIT)
#define ULONG_BITS_SQUARED (ULONG_BITS * ULONG_BITS)

typedef struct twod_square_graph {
  size_t my_nlocalverts_1d;
  size_t my_nlocalverts_2d_src, my_nlocalverts_2d_tgt;
  int64_t* nlocalverts_1d; /* Size is square_size**2 */
  int64_t* nlocalverts_2d; /* Size is square_size */
  int lg_nglobalverts;
  size_t nlocaledges;
  uint64_t* localedges;
  int my_row, my_col; /* For current rank */
  int square_side;
  int nvertexblocks; /* Global, in ULONG_BITS_SQUARED units */
  int this_rank_active; /* 0 if outside square */
  const tuple_graph* tg; /* Original graph used to build this one */
} twod_square_graph;

static inline int vertex_owner_2d(int64_t v, const twod_square_graph* g_ptr) {
  v /= ULONG_BITS_SQUARED;
  int64_t small_blocksize = g_ptr->nvertexblocks / g_ptr->square_side;
  int64_t large_blocksize = small_blocksize + 1;
  int block_size_cutoff = (int)(g_ptr->nvertexblocks % g_ptr->square_side);
  if (v < block_size_cutoff * large_blocksize) {
    return (int)(v / large_blocksize);
  } else {
    v -= block_size_cutoff * large_blocksize;
    return block_size_cutoff + (int)(v / small_blocksize);
  }
}

static inline int64_t vertex_offset_for_coord_2d(int coord, const twod_square_graph* g_ptr) {
  int64_t small_blocksize = g_ptr->nvertexblocks / g_ptr->square_side;
  int64_t large_blocksize = small_blocksize + 1;
  int block_size_cutoff = (int)(g_ptr->nvertexblocks % g_ptr->square_side);
  int64_t result = 0;
  if (coord < block_size_cutoff) {
    result = ((int64_t)coord * large_blocksize) * ULONG_BITS_SQUARED;
  } else {
    result = ((int64_t)block_size_cutoff * large_blocksize + (coord - block_size_cutoff) * small_blocksize) * ULONG_BITS_SQUARED;
  }
  return int64_min(result, (INT64_C(1) << g_ptr->lg_nglobalverts));
}

static inline size_t vertex_local_2d(int64_t v, const twod_square_graph* g_ptr) {
  int owner = vertex_owner_2d(v, g_ptr);
  v -= vertex_offset_for_coord_2d(owner, g_ptr);
  assert (v >= 0);
  return (size_t)v;
}

static inline int64_t vertex_to_global_2d(int coord, size_t vl, const twod_square_graph* g_ptr) {
  return (int64_t)vl + vertex_offset_for_coord_2d(coord, g_ptr);
}

static inline int vertex_owner_1d_from_2d_local(int owner_2d, size_t v, const twod_square_graph* g_ptr) {
  v /= ULONG_BITS_SQUARED;
  int64_t small_blocksize = g_ptr->nlocalverts_2d[owner_2d] / ULONG_BITS_SQUARED / g_ptr->square_side;
  int64_t large_blocksize = small_blocksize + 1;
  int block_size_cutoff = (int)((g_ptr->nlocalverts_2d[owner_2d] / ULONG_BITS_SQUARED) % g_ptr->square_side);
  if (v < block_size_cutoff * large_blocksize) {
    return (int)(v / large_blocksize);
  } else {
    v -= block_size_cutoff * large_blocksize;
    return block_size_cutoff + (int)(v / small_blocksize);
  }
}

static inline int vertex_owner_1d(int64_t v, const twod_square_graph* g_ptr) {
  int owner_2d = vertex_owner_2d(v, g_ptr);
  v -= vertex_offset_for_coord_2d(owner_2d, g_ptr);
  return owner_2d * g_ptr->square_side + vertex_owner_1d_from_2d_local(owner_2d, v, g_ptr);
}

static inline int64_t vertex_offset_for_coord_1d(int r, const twod_square_graph* g_ptr) {
  int row = r / g_ptr->square_side;
  int col = r % g_ptr->square_side;
  assert ((row < g_ptr->square_side && col < g_ptr->square_side) || (row == g_ptr->square_side && col == 0));
  if (row == g_ptr->square_side) {
    --row;
    col = g_ptr->square_side; /* Row needs to be within bounds, but col doesn't need to */
  }
  int64_t small_blocksize = g_ptr->nlocalverts_2d[row] / ULONG_BITS_SQUARED / g_ptr->square_side;
  int64_t large_blocksize = small_blocksize + 1;
  int block_size_cutoff = (int)((g_ptr->nlocalverts_2d[row] / ULONG_BITS_SQUARED) % g_ptr->square_side);
  int64_t row_part = vertex_offset_for_coord_2d(row, g_ptr);
  int64_t result = row_part;
  if (col < block_size_cutoff) {
    result += ((int64_t)col * large_blocksize) * ULONG_BITS_SQUARED;
  } else {
    result += ((int64_t)block_size_cutoff * large_blocksize + (col - block_size_cutoff) * small_blocksize) * ULONG_BITS_SQUARED;
  }
  return int64_min(result, (INT64_C(1) << g_ptr->lg_nglobalverts));
}

static inline size_t vertex_local_1d(int64_t v, const twod_square_graph* g_ptr) {
  return v - vertex_offset_for_coord_1d(vertex_owner_1d(v, g_ptr), g_ptr);
}

static inline int64_t vertex_to_global_1d(int rank, size_t vl, const twod_square_graph* g_ptr) {
  return vertex_offset_for_coord_1d(rank, g_ptr) + vl;
}

#define EDGE_OWNER(src, tgt, g_ptr) ((int)((vertex_owner_2d((src), (g_ptr)) * (g_ptr)->square_side) + vertex_owner_2d((tgt), (g_ptr))))
#define EDGE_LOCAL(src, tgt, g_ptr) ((uint64_t)(vertex_local_2d((src), (g_ptr)) << 32) + (uint64_t)(vertex_local_2d((tgt), (g_ptr))))

void convert_graph_to_twod_square(const tuple_graph* const tg, twod_square_graph* const g);
void free_twod_square_graph(twod_square_graph* const g);

#endif /* TWOD_SQUARE_H */
