/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#if !defined(RMAT_HEADER_)
#define RMAT_HEADER_

#include "generator/graph_generator.h"

/** Fill IJ with a randomly permuted R-MAT generated edge list. */
void rmat_edgelist (struct packed_edge *IJ, int64_t nedge, int SCALE,
		    double A, double B, double C);
void permute_vertex_labels (struct packed_edge * IJ, int64_t nedge, int64_t max_nvtx,
			    void * st, int64_t * newlabel);
void permute_edgelist (struct packed_edge * IJ, int64_t nedge, void *st);

#endif /* RMAT_HEADER_ */
