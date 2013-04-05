/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* See COPYING for license. */
#if !defined(GRAPH500_IMPL_HEADER_)
#define GRAPH500_IMPL_HEADER_

/** Pass the edge list to an external graph creation routine. */
int create_graph_from_edgelist (struct packed_edge *IJ, int64_t nedge,
				int64_t nv);

/** Create the BFS tree from a given source vertex. */
int make_bfs_tree (int64_t *bfs_tree_out, int64_t srcvtx);

/** Clean up. */
void destroy_graph (void);

#endif /* GRAPH500_IMPL_HEADER_ */
