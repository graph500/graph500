/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef ONED_CSC_H
#define ONED_CSC_H

#include "common.h"
#include <limits.h>

#define ULONG_BITS (sizeof(unsigned long) * CHAR_BIT)

typedef struct oned_csc_graph {
  size_t nlocalverts;
  int64_t max_nlocalverts;
  size_t nlocaledges;
  int lg_nglobalverts;
  int64_t nglobalverts;
  size_t *rowstarts;
  int64_t *column;
  int lg_local_queue_size; /* Queue size for replicated CSC BFS */
  const tuple_graph* tg; /* Original graph used to build this one */
} oned_csc_graph;

#define SWIZZLE_VERTEX(c) (((VERTEX_OWNER(c) << lg_local_queue_size) * ULONG_BITS) | VERTEX_LOCAL(c))
#define UNSWIZZLE_VERTEX(c) (MUL_SIZE((c) & ((INT64_C(1) << lg_local_queue_size) * ULONG_BITS - 1)) | ((((c) / ULONG_BITS) >> lg_local_queue_size)))

void convert_graph_to_oned_csc(const tuple_graph* const tg, oned_csc_graph* const g);
void free_oned_csc_graph(oned_csc_graph* const g);

#endif /* ONED_CSC_H */
