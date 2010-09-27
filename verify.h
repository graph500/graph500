/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#if !defined(VERIFY_HEADER_)
#define VERIFY_HEADER_

/** Verify a BFS tree, return volume or -1 if failed. */
int64_t verify_bfs_tree (int64_t *bfs_tree, int64_t max_bfsvtx,
			 int64_t root,
			 const int64_t *IJ, int64_t nedge);

#endif /* VERIFY_HEADER_ */
