/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include "common.h"
#include "oned_csr.h"
#include "redistribute.h"
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

typedef struct temp_csr_graph {
  size_t* restrict rowstarts;
  int64_t* restrict column;
  size_t nlocalverts;
  size_t nlocaledges;
  size_t nlocaledges_allocated; /* Actual size of column */
  int lg_nglobalverts;
} temp_csr_graph;

static void make_empty_csr(temp_csr_graph* restrict const outg /* All fields NULL or 0 */) {
  outg->rowstarts = (size_t*)xcalloc(1, sizeof(size_t));
  outg->column = NULL; /* Realloc can enlarge a NULL pointer */
  outg->nlocalverts = outg->nlocaledges = outg->nlocaledges_allocated = 0;
  outg->lg_nglobalverts = -1;
}

static void make_csr(const packed_edge* restrict const inbuf, temp_csr_graph* restrict const outg /* Must have memory and nlocalverts/nlocaledges filled in */) {
  size_t nrows = outg->nlocalverts;
  size_t inbuf_size = outg->nlocaledges;
  size_t* temp = (size_t*)xmalloc(nrows * sizeof(size_t));
  size_t* restrict rowstarts = outg->rowstarts;
  int64_t* restrict column = outg->column;
  {
    size_t* restrict counts = temp;
    memset(counts, 0, nrows * sizeof(size_t));
    ptrdiff_t i;
#pragma omp parallel for
    for (i = 0; i < (ptrdiff_t)inbuf_size; ++i) {
      assert ((size_t)(VERTEX_LOCAL(get_v0_from_edge(&inbuf[i]))) < nrows);
#pragma omp atomic
      ++counts[VERTEX_LOCAL(get_v0_from_edge(&inbuf[i]))];
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
      int64_t v1 = get_v1_from_edge(&inbuf[i]);
      assert ((size_t)(VERTEX_LOCAL(v0)) < nrows);
      size_t pos = __sync_fetch_and_add(&inserts[VERTEX_LOCAL(v0)], 1);
      assert (pos < inbuf_size);
      column[pos] = v1;
    }
  }
  free(temp);
}

/* Do merge: b = b union a */
static void merge_csr(temp_csr_graph* restrict const b,
                      const temp_csr_graph* restrict const a) {
  size_t a_nlocalverts = a->nlocalverts;
  size_t b_nlocalverts = b->nlocalverts;
  size_t a_nlocaledges = a->nlocaledges;
  size_t b_nlocaledges = b->nlocaledges;
  if (a->nlocalverts > b->nlocalverts) {
    ptrdiff_t old_b_nlocalverts = b_nlocalverts, i;
    b->rowstarts = (size_t*)xrealloc(b->rowstarts, (a_nlocalverts + 1) * sizeof(size_t));
    b_nlocalverts = b->nlocalverts = a->nlocalverts;
#pragma omp parallel for
    for (i = old_b_nlocalverts; i < b_nlocalverts; ++i) {
      b->rowstarts[i + 1] = b_nlocaledges;
    }
    b->lg_nglobalverts = a->lg_nglobalverts;
  }

  if (b_nlocaledges + a_nlocaledges > b->nlocaledges_allocated) {
    size_t new_alloc = b_nlocaledges + a_nlocaledges + (1 << 16);
    b->nlocaledges_allocated = new_alloc;
    b->column = (int64_t*)xrealloc(b->column, new_alloc * sizeof(int64_t));
  }
  memmove(&b->column[b->rowstarts[a_nlocalverts] + a_nlocaledges],
          &b->column[b->rowstarts[a_nlocalverts]],
          (b_nlocaledges - b->rowstarts[a_nlocalverts]) * sizeof(int64_t));
  ptrdiff_t i_plus_1;
  for (i_plus_1 = a_nlocalverts; i_plus_1 > 0; --i_plus_1) {
    ptrdiff_t i = i_plus_1 - 1;
    memmove(&b->column[b->rowstarts[i] + a->rowstarts[i]],
            &b->column[b->rowstarts[i]],
            (b->rowstarts[i + 1] - b->rowstarts[i]) * sizeof(int64_t));
    memcpy(&b->column[b->rowstarts[i + 1] + a->rowstarts[i]],
           &a->column[a->rowstarts[i]],
           (a->rowstarts[i + 1] - a->rowstarts[i]) * sizeof(int64_t));
  }
  b_nlocaledges = b->nlocaledges = b_nlocaledges + a_nlocaledges;
  ptrdiff_t i;
#pragma omp parallel for
  for (i = 0; i <= a_nlocalverts; ++i) {
    b->rowstarts[i] += a->rowstarts[i];
  }
#pragma omp parallel for if(a_nlocalverts != b_nlocalverts)
  for (i = a_nlocalverts + 1; i <= b_nlocalverts; ++i) {
    b->rowstarts[i] += a_nlocaledges;
  }
}

#define CONV1D_FUNCNAME \
  convert_graph_to_oned_csr_helper

#define CONV1D_EXTRA_PARAMS \
  oned_csr_graph* const g

#define CONV1D_DECLARE_AND_INIT_GRAPH_SO_FAR \
  temp_csr_graph graph_so_far = {NULL, NULL, 0, 0}; \
  make_empty_csr(&graph_so_far);

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
  temp_csr_graph t = { \
    /* rowstarts   */           (size_t*)xmalloc((size_t)(nlocalverts_so_far + 1) * sizeof(size_t)), \
    /* column      */           (int64_t*)xmalloc((size_t)(EDGES_RECEIVED_THIS_BLOCK) * sizeof(int64_t)), \
    /* nlocalverts */           (size_t)(nlocalverts_so_far), \
    /* nlocaledges */           (size_t)(EDGES_RECEIVED_THIS_BLOCK), \
    /* nlocaledges_allocated */ (size_t)(EDGES_RECEIVED_THIS_BLOCK), \
    /* lg_nglobalverts */       (int)(LG_NGLOBALVERTS_SO_FAR) \
  }; \
  make_csr((EDGES_TO_RECV), &t);

#define CONV1D_MERGE_INTO_GRAPH_SO_FAR \
  size_t new_alloc = graph_so_far.nlocaledges + edges_received_this_block * (block_count - ITERATE_TUPLE_GRAPH_BLOCK_NUMBER); \
  if (new_alloc > graph_so_far.nlocaledges_allocated) { \
    size_t new_alloc_real = new_alloc + (1 << 16); \
    graph_so_far.nlocaledges_allocated = new_alloc_real; \
    graph_so_far.column = (int64_t*)xrealloc(graph_so_far.column, new_alloc_real * sizeof(int64_t)); \
  } \
  merge_csr(&graph_so_far, &t);

#define CONV1D_FREE_PRECOMPRESSED_DATA \
  free(t.rowstarts); \
  free(t.column);

#define CONV1D_BUILD_FINAL_DATA_STRUCTURE_FROM_GRAPH_SO_FAR \
  g->nlocaledges = graph_so_far.nlocaledges; \
  g->rowstarts = graph_so_far.rowstarts; \
  g->column = (int64_t*)xrealloc(graph_so_far.column, (size_t)g->nlocaledges * sizeof(int64_t)); \
  size_t nlocalverts = graph_so_far.nlocalverts; \
  g->nlocalverts = nlocalverts; \
  g->max_nlocalverts = nlocalverts; /* Now same on all ranks */ \
  g->lg_nglobalverts = graph_so_far.lg_nglobalverts; \
  g->nglobalverts = INT64_C(1) << graph_so_far.lg_nglobalverts;

#define CONV1D_CLEAR_GRAPH_SO_FAR \
  free(graph_so_far.rowstarts); graph_so_far.rowstarts = NULL; \
  free(graph_so_far.column); graph_so_far.column = NULL; \
  graph_so_far.nlocalverts = graph_so_far.nlocaledges = graph_so_far.nlocaledges_allocated = 0;
  
static MAKE_REDISTRIBUTE_FUNC(CONV1D_FUNCNAME, CONV1D_EXTRA_PARAMS, CONV1D_DECLARE_AND_INIT_GRAPH_SO_FAR, CONV1D_CALL_ON_EDGES, CONV1D_EDGE_BUFFER_TYPE, CONV1D_EDGE_BUFFER_MPI_TYPE, CONV1D_PRECOMPRESS_INCOMING_DATA, CONV1D_MERGE_INTO_GRAPH_SO_FAR, CONV1D_FREE_PRECOMPRESSED_DATA, CONV1D_BUILD_FINAL_DATA_STRUCTURE_FROM_GRAPH_SO_FAR, CONV1D_CLEAR_GRAPH_SO_FAR)

void convert_graph_to_oned_csr(const tuple_graph* const tg, oned_csr_graph* const g) { \
  g->tg = tg;
  g->nlocaledges = 0;
  convert_graph_to_oned_csr_helper(tg, g);
}

void free_oned_csr_graph(oned_csr_graph* const g) {
  if (g->rowstarts != NULL) {free(g->rowstarts); g->rowstarts = NULL;}
  if (g->column != NULL) {free(g->column); g->column = NULL;}
}

