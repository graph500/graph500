/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef ONED_CSR_H
#define ONED_CSR_H

#include "common.h"

typedef struct oned_csr_graph {
  size_t nlocalverts;
  int64_t max_nlocalverts;
  size_t nlocaledges;
  int lg_nglobalverts;
  int64_t nglobalverts;
  size_t *rowstarts;
  int64_t *column;
  const tuple_graph* tg; /* Original graph used to build this one */
} oned_csr_graph;

void convert_graph_to_oned_csr(const tuple_graph* const tg, oned_csr_graph* const g);
void free_oned_csr_graph(oned_csr_graph* const g);

#endif /* ONED_CSR_H */
