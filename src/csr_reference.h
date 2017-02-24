/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine, Anton Korzh                                 */

#ifndef CSR_REFERENCE_H
#define CSR_REFERENCE_H

#include "common.h"

typedef struct oned_csr_graph {
	size_t nlocalverts;
	int64_t max_nlocalverts;
	size_t nlocaledges;
	int lg_nglobalverts;
	int64_t nglobalverts,notisolated;
	int *rowstarts;
	int64_t *column;
#ifdef SSSP 
	float *weights;
#endif
	const tuple_graph* tg;
} oned_csr_graph;

void convert_graph_to_oned_csr(const tuple_graph* const tg, oned_csr_graph* const g);
void free_oned_csr_graph(oned_csr_graph* const g);

//#define BYTES_PER_VERTEX 8
//#define COLUMN(i) column[i]
//#define SETCOLUMN(a,b) column[a]=b;
#define BYTES_PER_VERTEX 6
#define SETCOLUMN(a,b) memcpy(((char*)column)+(BYTES_PER_VERTEX*a),&b,BYTES_PER_VERTEX)
#define COLUMN(i) (*(int64_t*)(((char*)column)+(BYTES_PER_VERTEX*i)) & (int64_t)(0xffffffffffffffffULL>>(64-8*BYTES_PER_VERTEX)))

#endif /* CSR_REFERENCE_H */
