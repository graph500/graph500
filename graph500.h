/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#if !defined(GRAPH500_HEADER_)
#define GRAPH500_HEADER_

#define NAME "Graph500 sequential list"
#define VERSION 0

/** Pass the edge list to an external graph creation routine. */
int create_graph_from_edgelist (int64_t *IJ, int64_t nedge);

/** Create the BFS tree from a given source vertex. */
int make_bfs_tree (int64_t *bfs_tree_out, int64_t *max_vtx_out,
		   int64_t srcvtx);

/** Clean up. */
void destroy_graph (void);

#endif /* GRAPH500_HEADER_ */
