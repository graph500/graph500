/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#if !defined(GRAPH500_HEADER_)
#define GRAPH500_HEADER_

#define NAME "Graph500 sequential list"
#define VERSION 0

#include "generator/graph_generator.h"

/** Pass the edge list to an external graph creation routine. */
int create_graph_from_edgelist (struct packed_edge *IJ, int64_t nedge);

/** Create the BFS tree from a given source vertex. */
int make_bfs_tree (int64_t *bfs_tree_out, int64_t *max_vtx_out,
		   int64_t srcvtx);

/** Clean up. */
void destroy_graph (void);

#endif /* GRAPH500_HEADER_ */
