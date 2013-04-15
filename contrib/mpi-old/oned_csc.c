/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include "common.h"
#include "oned_csc.h"
#include "redistribute.h"
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

typedef struct temp_csc_graph {
  size_t* restrict rowstarts;
  int64_t* restrict column;
  size_t nlocalverts;
  int lg_nglobalverts;
  int64_t nglobalverts;
  size_t nlocaledges;
  size_t nlocaledges_allocated; /* Actual size of column */
  int lg_local_queue_size;
  size_t nrows; /* One less than size of rowstarts */
} temp_csc_graph;

static void make_empty_csc(temp_csc_graph* restrict const outg /* All fields NULL or 0 */) {
  outg->rowstarts = (size_t*)xcalloc(1, sizeof(size_t));
  outg->column = NULL; /* Realloc can enlarge a NULL pointer */
  outg->nlocalverts = outg->nglobalverts = outg->nlocaledges = outg->nlocaledges_allocated = 0;
  outg->lg_nglobalverts = -1;
  outg->lg_local_queue_size = -1;
  outg->nrows = 0;
}

static void make_csc(const packed_edge* restrict const inbuf, temp_csc_graph* restrict const outg /* Must have memory and nlocalverts/nglobalverts/nlocaledges filled in */) {
  size_t nrows = outg->nrows;
  size_t inbuf_size = outg->nlocaledges;
  size_t* temp = (size_t*)xmalloc(nrows * sizeof(size_t));
  size_t* restrict rowstarts = outg->rowstarts;
  int64_t* restrict column = outg->column;
  int lg_local_queue_size = outg->lg_local_queue_size;
  {
    size_t* restrict counts = temp;
    memset(counts, 0, nrows * sizeof(size_t));
    ptrdiff_t i;
#pragma omp parallel for
    for (i = 0; i < (ptrdiff_t)inbuf_size; ++i) {
      assert ((size_t)(SWIZZLE_VERTEX(get_v1_from_edge(&inbuf[i])) / ULONG_BITS) < nrows);
#pragma omp atomic
      ++counts[SWIZZLE_VERTEX(get_v1_from_edge(&inbuf[i])) / ULONG_BITS];
    }
    rowstarts[0] = 0;
    for (i = 0; i < nrows; ++i) {
      rowstarts[i + 1] = rowstarts[i] + counts[i];
    }
  }
  {
    size_t* restrict inserts = temp;
    memcpy(inserts, rowstarts, nrows * sizeof(size_t));
    ptrdiff_t i;
#pragma omp parallel for
    for (i = 0; i < (ptrdiff_t)inbuf_size; ++i) {
      int64_t v0 = get_v0_from_edge(&inbuf[i]);
      int64_t v1 = SWIZZLE_VERTEX(get_v1_from_edge(&inbuf[i]));
      // fprintf(stderr, "%d: Raw edge is (%" PRId64 ", %" PRId64 ") -> (%zu, %" PRId64 " = %" PRId64 ")\n", rank, v0, get_v1_from_edge(&inbuf[i]), VERTEX_LOCAL(v0), v1, UNSWIZZLE_VERTEX(v1));
      size_t pos = __sync_fetch_and_add(&inserts[(v1) / ULONG_BITS], 1);
      column[pos] = (v1 % ULONG_BITS) + VERTEX_LOCAL(v0) * ULONG_BITS;
      // fprintf(stderr, "%d: Stored as (row %" PRId64 ", col %" PRId64 "/%" PRId64 ")\n", rank, (v1) / ULONG_BITS, column[pos] % ULONG_BITS, column[pos] / ULONG_BITS);
    }
  }
  free(temp); temp = NULL;
}

/* Do merge: b = b union a */
static void merge_csc(temp_csc_graph* restrict const b,
                      temp_csc_graph* restrict const a) {
  if (b->lg_local_queue_size == -1) { // b is empty
    if (b->rowstarts != NULL) {free(b->rowstarts); b->rowstarts = NULL;}
    if (b->column != NULL) {free(b->column); b->column = NULL;}
    *b = *a;
    a->rowstarts = NULL;
    a->column = NULL;
    return;
  } else if (a->nglobalverts != b->nglobalverts) {
    /* Redistribution wrapper should restart in this case, not try to do a merge. */
    fprintf(stderr, "%d: a->nglobalverts=%" PRId64 " != b->nglobalverts=%" PRId64 "\n", rank, a->nglobalverts, b->nglobalverts);
    MPI_Abort(MPI_COMM_WORLD, 5);
  } else {
    assert (a->lg_local_queue_size == b->lg_local_queue_size);
    assert (a->nrows == b->nrows);
    assert (a->lg_nglobalverts == b->lg_nglobalverts);
    size_t a_nlocaledges = a->nlocaledges;
    size_t b_nlocaledges = b->nlocaledges;
    size_t nrows = b->nrows;

    if (b_nlocaledges + a_nlocaledges > b->nlocaledges_allocated) {
      size_t new_alloc = b_nlocaledges + a_nlocaledges + (1 << 16);
      b->nlocaledges_allocated = new_alloc;
      b->column = (int64_t*)xrealloc(b->column, new_alloc * sizeof(int64_t));
    }
    ptrdiff_t i_plus_1;
    /* This loop needs to be sequential. */
    for (i_plus_1 = nrows; i_plus_1 > 0; --i_plus_1) {
      ptrdiff_t i = i_plus_1 - 1;
      memmove(&b->column[b->rowstarts[i] + a->rowstarts[i]],
              &b->column[b->rowstarts[i]],
              (b->rowstarts[i + 1] - b->rowstarts[i]) * sizeof(int64_t));
    }
    /* This loop can be parallel. */
#pragma omp parallel for
    for (i_plus_1 = nrows; i_plus_1 > 0; --i_plus_1) {
      ptrdiff_t i = i_plus_1 - 1;
      memcpy(&b->column[b->rowstarts[i + 1] + a->rowstarts[i]],
             &a->column[a->rowstarts[i]],
             (a->rowstarts[i + 1] - a->rowstarts[i]) * sizeof(int64_t));
    }
    b_nlocaledges = b->nlocaledges = b_nlocaledges + a_nlocaledges;
    ptrdiff_t i;
#pragma omp parallel for
    for (i = 0; i <= nrows; ++i) {
      b->rowstarts[i] += a->rowstarts[i];
    }
    free(a->column); a->column = NULL;
    free(a->rowstarts); a->rowstarts = NULL;
  }
}

#define CONV1D_FUNCNAME \
  convert_graph_to_oned_csc_helper

#define CONV1D_EXTRA_PARAMS \
  oned_csc_graph* const g

#define CONV1D_DECLARE_AND_INIT_GRAPH_SO_FAR \
  temp_csc_graph graph_so_far = {NULL, NULL, 0, 0, 0}; \
  make_empty_csc(&graph_so_far);

#define CONV1D_CALL_ON_EDGES(V0, V1, LG_NGLOBALVERTS_SO_FAR, CONT) \
  CONT(VERTEX_OWNER((V0)), CONV1D_WRITE_EDGE_NORMAL) \
  CONT(VERTEX_OWNER((V1)), CONV1D_WRITE_EDGE_FLIPPED)
#define CONV1D_WRITE_EDGE_NORMAL(BUF, V0, V1) \
  write_edge(BUF, V0, V1);
#define CONV1D_WRITE_EDGE_FLIPPED(BUF, V0, V1) \
  write_edge(BUF, V1, V0);

#define CONV1D_EDGE_BUFFER_TYPE \
  packed_edge

#define CONV1D_EDGE_BUFFER_MPI_TYPE \
  packed_edge_mpi_type

#define CONV1D_PRECOMPRESS_INCOMING_DATA(LG_NGLOBALVERTS_SO_FAR, EDGES_TO_RECV, EDGES_RECEIVED_THIS_BLOCK) \
  size_t nlocalverts_so_far = (size_t)DIV_SIZE((UINT64_C(1) << (LG_NGLOBALVERTS_SO_FAR)) + size - 1); \
  size_t t_nrows = (size_t)(MUL_SIZE((nlocalverts_so_far + ULONG_BITS * ULONG_BITS - 1) / ULONG_BITS / ULONG_BITS * ULONG_BITS)); \
  temp_csc_graph t = { \
    /* rowstarts   */           (size_t*)xmalloc((t_nrows + 1) * sizeof(size_t)), \
    /* column      */           (int64_t*)xmalloc((size_t)(EDGES_RECEIVED_THIS_BLOCK) * sizeof(int64_t)), \
    /* nlocalverts */           (size_t)(nlocalverts_so_far), \
    /* lg_nglobalverts */       (int)(LG_NGLOBALVERTS_SO_FAR), \
    /* nglobalverts */          (int64_t)(INT64_C(1) << (LG_NGLOBALVERTS_SO_FAR)), \
    /* nlocaledges */           (size_t)(EDGES_RECEIVED_THIS_BLOCK), \
    /* nlocaledges_allocated */ (size_t)(EDGES_RECEIVED_THIS_BLOCK), \
    /* lg_local_queue_size */   -1, /* Filled in later */ \
    /* nrows */                 t_nrows \
  }; \
  { \
    t.lg_local_queue_size = lg_int64_t(DIV_SIZE(t_nrows)); \
    make_csc((EDGES_TO_RECV), &t); \
  }

#define CONV1D_MERGE_INTO_GRAPH_SO_FAR \
  size_t new_alloc = graph_so_far.nlocaledges + edges_received_this_block * (block_count - ITERATE_TUPLE_GRAPH_BLOCK_NUMBER); \
  if (graph_so_far.lg_local_queue_size != -1 && new_alloc > graph_so_far.nlocaledges_allocated) { \
    size_t new_alloc_real = new_alloc + (1 << 16); \
    graph_so_far.nlocaledges_allocated = new_alloc_real; \
    graph_so_far.column = (int64_t*)xrealloc(graph_so_far.column, new_alloc_real * sizeof(int64_t)); \
  } \
  merge_csc(&graph_so_far, &t);

#define CONV1D_FREE_PRECOMPRESSED_DATA \
  if (t.rowstarts != NULL) {free(t.rowstarts); t.rowstarts = NULL;} \
  if (t.column != NULL) {free(t.column); t.column = NULL;}

#define CONV1D_BUILD_FINAL_DATA_STRUCTURE_FROM_GRAPH_SO_FAR \
  g->nlocaledges = graph_so_far.nlocaledges; \
  g->rowstarts = graph_so_far.rowstarts; \
  graph_so_far.rowstarts = NULL; \
  g->column = (int64_t*)xrealloc(graph_so_far.column, (size_t)g->nlocaledges * sizeof(int64_t)); \
  graph_so_far.column = NULL; \
  g->lg_local_queue_size = graph_so_far.lg_local_queue_size; \
  size_t nlocalverts = graph_so_far.nlocalverts; \
  g->nlocalverts = nlocalverts; \
  g->max_nlocalverts = nlocalverts; /* Now same on all ranks */ \
  g->lg_nglobalverts = graph_so_far.lg_nglobalverts; \
  g->nglobalverts = INT64_C(1) << graph_so_far.lg_nglobalverts;

#define CONV1D_CLEAR_GRAPH_SO_FAR \
  free(graph_so_far.rowstarts); graph_so_far.rowstarts = NULL; \
  free(graph_so_far.column); graph_so_far.column = NULL; \
  graph_so_far.nlocalverts = graph_so_far.nlocaledges = graph_so_far.nlocaledges_allocated = 0; \
  graph_so_far.lg_local_queue_size = -1; \
  graph_so_far.nrows = 0;
  
static MAKE_REDISTRIBUTE_FUNC(CONV1D_FUNCNAME, CONV1D_EXTRA_PARAMS, CONV1D_DECLARE_AND_INIT_GRAPH_SO_FAR, CONV1D_CALL_ON_EDGES, CONV1D_EDGE_BUFFER_TYPE, CONV1D_EDGE_BUFFER_MPI_TYPE, CONV1D_PRECOMPRESS_INCOMING_DATA, CONV1D_MERGE_INTO_GRAPH_SO_FAR, CONV1D_FREE_PRECOMPRESSED_DATA, CONV1D_BUILD_FINAL_DATA_STRUCTURE_FROM_GRAPH_SO_FAR, CONV1D_CLEAR_GRAPH_SO_FAR)

void convert_graph_to_oned_csc(const tuple_graph* const tg, oned_csc_graph* const g) { \
  g->tg = tg;
  g->nlocaledges = 0;
  convert_graph_to_oned_csc_helper(tg, g);
  g->max_nlocalverts = (int64_t)(g->nlocalverts);
  MPI_Allreduce(MPI_IN_PLACE, &g->max_nlocalverts, 1, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);
  int64_t local_queue_summary_size = (g->max_nlocalverts + ULONG_BITS * ULONG_BITS - 1) / ULONG_BITS / ULONG_BITS;
  int64_t local_queue_size = local_queue_summary_size * ULONG_BITS;
  if (g->lg_local_queue_size != lg_int64_t(local_queue_size)) {
    fprintf(stderr, "%d: lg_local_queue_size mismatch: graph redistribution computed %d, convert_graph_to_oned_csc outer computed %d from %" PRId64 "\n", rank, g->lg_local_queue_size, lg_int64_t(local_queue_size), local_queue_size);
    MPI_Abort(MPI_COMM_WORLD, 6);
  }
}

void free_oned_csc_graph(oned_csc_graph* const g) {
  if (g->rowstarts != NULL) {free(g->rowstarts); g->rowstarts = NULL;}
  if (g->column != NULL) {free(g->column); g->column = NULL;}
}

