/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _FILE_OFFSET_BITS 64
#define _THREAD_SAFE
#define _XOPEN_SOURCE 600
#define _XOPEN_SOURCE_EXTENDED
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <assert.h>

#include "compat.h"
#include "xalloc.h"

static int
compute_levels (int64_t * level,
		int64_t nv, const int64_t * restrict bfs_tree, int64_t root)
{
  int err = 0;
  int64_t k;
  for (k = 0; k < nv; ++k) level[k] = -1;
  level[root] = 0;

  for (k = 0; k < nv; ++k) {
    if (level[k] >= 0) continue;
    if (!err && bfs_tree[k] >= 0 && k != root) {
      int64_t parent = k;
      int64_t nhop = 0;
      /* Run up the tree until we encounter an already-leveled vertex. */
      while (parent >= 0 && level[parent] < 0 && nhop < nv) {
	assert (parent != bfs_tree[parent]);
	parent = bfs_tree[parent];
	++nhop;
      }
      if (nhop >= nv) err = -1; /* Cycle. */
      if (parent < 0) err = -2; /* Ran off the end. */

      /* Now assign levels until we meet an already-leveled vertex */
      /* NOTE: This permits benign races if parallelized. */
      nhop += level[parent];
      parent = k;
      while (level[parent] < 0) {
	assert (nhop > 0);
	level[parent] = nhop--;
	parent = bfs_tree[parent];
      }
      assert (nhop == level[parent]);

      /* Internal check to catch mistakes in races... */
#if !defined(NDEBUG)
      nhop = 0;
      parent = k;
      int64_t lastlvl = level[k]+1;
      while (level[parent] > 0) {
	assert (lastlvl == 1 + level[parent]);
	lastlvl = level[parent];
	parent = bfs_tree[parent];
	++nhop;
      }
#endif
    }
  }
  return err;
}

int64_t
verify_bfs_tree (int64_t *bfs_tree_in, int64_t max_bfsvtx,
		 int64_t root,
		 const int64_t *IJ_in, int64_t nedge)
{
  int64_t * restrict bfs_tree = bfs_tree_in;
  const int64_t * restrict IJ = IJ_in;

  int err, nedge_traversed;
  int64_t * restrict seen_edge, * restrict level;

  const int64_t nv = max_bfsvtx+1;
  int64_t k;

  /*
    This code is horrifically contorted because many compilers
    complain about continue, return, etc. in parallel sections.
  */

  if (root > max_bfsvtx || bfs_tree[root] != root)
    return -999;

  err = 0;
  nedge_traversed = 0;
  seen_edge = xmalloc_large (2 * (nv) * sizeof (*seen_edge));
  level = &seen_edge[nv];

  for (k = 0; k < nv; ++k)
    seen_edge[k] = 0;

  err = compute_levels (level, nv, bfs_tree, root);

  if (err) goto done;

  for (k = 0; k < 2*nedge; k+=2) {
    const int64_t i = IJ[k];
    const int64_t j = IJ[k+1];
    int64_t lvldiff;

    if (i > max_bfsvtx && j <= max_bfsvtx) err = -10;
    if (j > max_bfsvtx && i <= max_bfsvtx) err = -11;
    if (err || i > max_bfsvtx /* both i & j are on the same side of max_bfsvtx */)
      continue;

    /* All neighbors must be in the tree. */
    if (bfs_tree[i] >= 0 && bfs_tree[j] < 0) err = -12;
    if (bfs_tree[j] >= 0 && bfs_tree[i] < 0) err = -13;
    if (err || bfs_tree[i] < 0 /* both i & j have the same sign */)
      continue;

    /* Both i and j are in the tree, count as a traversed edge.

       NOTE: This counts self-edges and repeated edges.  They're
       part of the input data.
    */
    ++nedge_traversed;
    /* Mark seen tree edges. */
    if (i != j) {
      if (bfs_tree[i] == j)
	seen_edge[i] = 1;
      if (bfs_tree[j] == i)
	seen_edge[j] = 1;
    }
    lvldiff = level[i] - level[j];
    /* Check that the levels differ by no more than one. */
    if (lvldiff > 1 || lvldiff < -1)
      err = -14;
  }

  if (err) goto done;

  /* Check that every BFS edge was seen and that there's only one root. */
  for (k = 0; k < nv; ++k)
    if (k != root) {
      if (bfs_tree[k] >= 0 && !seen_edge[k])
	err = -15;
      if (bfs_tree[k] == k)
	err = -16;
    }

 done:

  xfree_large (seen_edge);
  if (err) return err;
  return nedge_traversed;
}
