/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include "common.h"
#include "twod_square.h"
#include "redistribute.h"
#include <math.h>
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

typedef struct temp_graph_storage {
  uint64_t* edges;
  size_t nlocaledges;
  size_t nlocaledges_allocated;
  int lg_nglobalverts;
} temp_graph_storage;

static void make_temp_graph(const uint64_t* restrict const inbuf, temp_graph_storage* restrict const outg /* Must have memory and nlocaledges filled in */) {
  assert (outg->nlocaledges <= outg->nlocaledges_allocated);
  memcpy(outg->edges, inbuf, outg->nlocaledges * sizeof(uint64_t));
}

/* Do merge: b = b union a */
static void merge_temp_graph(temp_graph_storage* restrict const b,
                             const temp_graph_storage* restrict const a) {
  if (b->nlocaledges_allocated < a->nlocaledges + b->nlocaledges) {
    b->nlocaledges_allocated = a->nlocaledges + b->nlocaledges;
    b->edges = (uint64_t*)xrealloc(b->edges, b->nlocaledges_allocated);
  }
  memcpy(b->edges + b->nlocaledges, a->edges, a->nlocaledges * sizeof(uint64_t));
  b->nlocaledges += a->nlocaledges;
  if (a->lg_nglobalverts > b->lg_nglobalverts) b->lg_nglobalverts = a->lg_nglobalverts;
}

#define CONV2D_FUNCNAME \
  convert_graph_to_twod_square_helper

#define CONV2D_EXTRA_PARAMS \
  twod_square_graph* const g

#define CONV2D_DECLARE_AND_INIT_GRAPH_SO_FAR \
  temp_graph_storage graph_so_far = {NULL, 0, 0, -1};

#define CONV2D_CALL_ON_EDGES(V0, V1, LG_NGLOBALVERTS_SO_FAR, CONT) \
  (g)->lg_nglobalverts = (LG_NGLOBALVERTS_SO_FAR); \
  size_t nvertexblocks = ((INT64_C(1) << (LG_NGLOBALVERTS_SO_FAR)) + ULONG_BITS_SQUARED - 1) / ULONG_BITS_SQUARED; \
  (g)->nvertexblocks = nvertexblocks; \
  CONT(EDGE_OWNER((V0), (V1), (g)), CONV2D_WRITE_EDGE_NORMAL) \
  CONT(EDGE_OWNER((V1), (V0), (g)), CONV2D_WRITE_EDGE_FLIPPED)
#define CONV2D_WRITE_EDGE_NORMAL(BUF, V0, V1) \
  *(BUF) = EDGE_LOCAL((V0), (V1), (g));
#define CONV2D_WRITE_EDGE_FLIPPED(BUF, V0, V1) \
  *(BUF) = EDGE_LOCAL((V1), (V0), (g));

#define CONV2D_EDGE_BUFFER_TYPE \
  uint64_t

#define CONV2D_EDGE_BUFFER_MPI_TYPE \
  MPI_UINT64_T

#define CONV2D_PRECOMPRESS_INCOMING_DATA(LG_NGLOBALVERTS_SO_FAR, EDGES_TO_RECV, EDGES_RECEIVED_THIS_BLOCK) \
  temp_graph_storage t = { \
    /* edges                 */ (uint64_t*)xmalloc((size_t)(EDGES_RECEIVED_THIS_BLOCK) * sizeof(uint64_t)), \
    /* nlocaledges           */ (size_t)(EDGES_RECEIVED_THIS_BLOCK), \
    /* nlocaledges_allocated */ (size_t)(EDGES_RECEIVED_THIS_BLOCK), \
    /* lg_nglobalverts       */ (int)(LG_NGLOBALVERTS_SO_FAR) \
  }; \
  make_temp_graph((EDGES_TO_RECV), &t);

#define CONV2D_MERGE_INTO_GRAPH_SO_FAR \
  size_t new_alloc = graph_so_far.nlocaledges + edges_received_this_block * (block_count - ITERATE_TUPLE_GRAPH_BLOCK_NUMBER); \
  if (new_alloc > graph_so_far.nlocaledges_allocated) { \
    size_t new_alloc_real = new_alloc + (1 << 16); \
    graph_so_far.nlocaledges_allocated = new_alloc_real; \
    graph_so_far.edges = (uint64_t*)xrealloc(graph_so_far.edges, new_alloc_real * sizeof(uint64_t)); \
  } \
  merge_temp_graph(&graph_so_far, &t);

#define CONV2D_FREE_PRECOMPRESSED_DATA \
  free(t.edges);

static inline void build_final_data_structure(twod_square_graph* g, temp_graph_storage* graph_so_far) {
  if (!g->this_rank_active) {g->localedges = NULL; g->nlocalverts_1d = NULL; g->nlocalverts_2d = NULL; return;}
  int64_t nglobalverts = INT64_C(1) << graph_so_far->lg_nglobalverts;
  size_t nvertexblocks = (nglobalverts + ULONG_BITS_SQUARED - 1) / ULONG_BITS_SQUARED;
  g->nlocalverts_2d = (int64_t*)xmalloc(g->square_side * sizeof(int64_t));
  {
    int i;
    for (i = 0; i < g->square_side; ++i) {
      g->nlocalverts_2d[i] =
        int64_min(
          (nvertexblocks / g->square_side + (i < nvertexblocks % g->square_side)) * (int64_t)ULONG_BITS_SQUARED,
          nglobalverts);
    }
  }
  g->nlocalverts_1d = (int64_t*)xmalloc(g->square_side * g->square_side * sizeof(int64_t));
  {
    int i;
    for (i = 0; i < g->square_side * g->square_side; ++i) {
      size_t nvertexblocks_this_row = nvertexblocks / g->square_side + ((i / g->square_side) < nvertexblocks % g->square_side);
      g->nlocalverts_1d[i] =
        int64_min(
          (nvertexblocks_this_row / g->square_side + ((i % g->square_side) < nvertexblocks_this_row % g->square_side)) * (int64_t)ULONG_BITS_SQUARED,
          nglobalverts);
    }
  }
  size_t nmyvertexblocks_1d = g->nlocalverts_1d[rank] / ULONG_BITS_SQUARED;
  g->localedges = (uint64_t*)xrealloc(graph_so_far->edges, graph_so_far->nlocaledges * sizeof(uint64_t));
  graph_so_far->edges = NULL;
  g->nlocaledges = graph_so_far->nlocaledges;
  g->lg_nglobalverts = graph_so_far->lg_nglobalverts;
  g->my_nlocalverts_1d = int64_min(nmyvertexblocks_1d * ULONG_BITS_SQUARED, INT64_C(1) << graph_so_far->lg_nglobalverts);
  g->my_nlocalverts_2d_src = g->nlocalverts_2d[g->my_row];
  size_t nmyvertexblocks_2d_row = nvertexblocks / g->square_side + (g->my_row < nvertexblocks % g->square_side);
  assert (g->my_nlocalverts_2d_src == nmyvertexblocks_2d_row * ULONG_BITS_SQUARED);
  g->my_nlocalverts_2d_tgt = g->nlocalverts_2d[g->my_col];
  g->nvertexblocks = nvertexblocks;
  assert (vertex_offset_for_coord_1d(rank + 1, g) - vertex_offset_for_coord_1d(rank, g) == g->my_nlocalverts_1d);
}

#define CONV2D_BUILD_FINAL_DATA_STRUCTURE_FROM_GRAPH_SO_FAR \
  build_final_data_structure(g, &graph_so_far);

#define CONV2D_CLEAR_GRAPH_SO_FAR \
  graph_so_far.edges = (uint64_t*)xrealloc(graph_so_far.edges, 0); \
  graph_so_far.nlocaledges = graph_so_far.nlocaledges_allocated = 0;
  
static MAKE_REDISTRIBUTE_FUNC(CONV2D_FUNCNAME, CONV2D_EXTRA_PARAMS, CONV2D_DECLARE_AND_INIT_GRAPH_SO_FAR, CONV2D_CALL_ON_EDGES, CONV2D_EDGE_BUFFER_TYPE, CONV2D_EDGE_BUFFER_MPI_TYPE, CONV2D_PRECOMPRESS_INCOMING_DATA, CONV2D_MERGE_INTO_GRAPH_SO_FAR, CONV2D_FREE_PRECOMPRESSED_DATA, CONV2D_BUILD_FINAL_DATA_STRUCTURE_FROM_GRAPH_SO_FAR, CONV2D_CLEAR_GRAPH_SO_FAR)

void convert_graph_to_twod_square(const tuple_graph* const tg, twod_square_graph* const g) { \
  g->tg = tg;
  g->nlocaledges = 0;
  g->square_side = (int)floor(sqrt(size));
  if (rank >= g->square_side * g->square_side) {
    g->this_rank_active = 0;
  } else {
    g->this_rank_active = 1;
    g->my_row = rank / g->square_side;
    g->my_col = rank % g->square_side;
  }
  convert_graph_to_twod_square_helper(tg, g);
}

void free_twod_square_graph(twod_square_graph* const g) {
  if (g->localedges != NULL) {free(g->localedges); g->localedges = NULL;}
  if (g->nlocalverts_1d != NULL) {free(g->nlocalverts_1d); g->nlocalverts_1d = NULL;}
  if (g->nlocalverts_2d != NULL) {free(g->nlocalverts_2d); g->nlocalverts_2d = NULL;}
}

