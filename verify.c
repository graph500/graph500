/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010-2011,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <assert.h>

#include "globals.h"
#include "xalloc.h"
#include "packed_edge.h"
#include "generator.h"
#include "verify.h"

static int
compute_levels (int64_t * level,
		int64_t nv, const int64_t * restrict bfs_tree, int64_t root)
{
  int err = 0;

  OMP("omp parallel shared(err)") {
    int terr;
    int64_t k;

    OMP("omp for")
      for (k = 0; k < nv; ++k)
	level[k] = (k == root? 0 : -1);

    OMP("omp for") MTA("mta assert parallel") MTA("mta use 100 streams")
      for (k = 0; k < nv; ++k) {
	if (level[k] >= 0) continue;
	terr = err;
	if (!terr && bfs_tree[k] >= 0 && k != root) {
	  int64_t parent = k;
	  int64_t nhop = 0;
	  /* Run up the tree until we encounter an already-leveled vertex. */
	  while (parent >= 0 && level[parent] < 0 && nhop < nv) {
	    assert (parent != bfs_tree[parent]);
	    parent = bfs_tree[parent];
	    ++nhop;
	  }
	  if (nhop >= nv) terr = -1; /* Cycle. */
	  if (parent < 0) terr = -2; /* Ran off the end. */

	  if (!terr) {
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
	if (terr) { err = terr;	OMP("omp flush (err)"); }
      }
  }
  return err;
}

int64_t
verify_bfs_tree (const int64_t *bfs_tree_in, const int64_t *level_in,
		 const int64_t root,
		 const struct packed_edge *IJ_in, int64_t nedge)
{
  const int64_t * restrict bfs_tree = bfs_tree_in;
  const int64_t * restrict level = level_in;
#if defined(STORED_EDGELIST)
  const struct packed_edge * restrict IJ = IJ_in;
#endif

  int err;
  int64_t maxdepth = -1;
  int64_t * restrict seen_edge;

  /*
    This code is horrifically contorted because many compilers
    complain about continue, return, etc. in parallel sections.
  */

  if (root >= NV || bfs_tree[root] != root)
    return -999;

  err = 0;
  if (!level) {
    int64_t * restrict level_tmp;
    seen_edge = xmalloc_large (2 * NV * sizeof (*seen_edge));
    level_tmp = &seen_edge[NV];
    err = compute_levels (level_tmp, NV, bfs_tree, root);
    level = level_tmp; /* Type launder to const. */
  } else {
    seen_edge = xmalloc_large (NV * sizeof (*seen_edge));
  }

  if (err) goto done;

  OMP("omp parallel shared(err)") {
    int64_t k;
    int terr = 0;
    int64_t tmaxdepth = -1;
    OMP("omp for")
      for (k = 0; k < NV; ++k)
	seen_edge[k] = 0;

    OMP("omp for")
    MTA("mta assert parallel") MTA("mta use 100 streams")
      for (k = 0; k < nedge; ++k) {
	int64_t i, j;
	uint8_t w;
#if !defined(STORED_EDGELIST)
	make_edge (k, &i, &j, &w);
#else
	i = get_v0_from_edge(&IJ[k]);
	j = get_v1_from_edge(&IJ[k]);
#endif
	int64_t lvldiff;
	terr = err;

	if (i < 0 || j < 0) continue;
	if (i >= NV && j < NV) terr = -10;
	if (j >= NV && i < NV) terr = -11;
	if (terr) { err = terr; OMP("omp flush(err)"); }
	if (terr || i >= NV /* both i & j are on the same side of NV */)
	  continue;

	/* All neighbors must be in the tree. */
	if (bfs_tree[i] >= 0 && bfs_tree[j] < 0) terr = -12;
	if (bfs_tree[j] >= 0 && bfs_tree[i] < 0) terr = -13;
	if (terr) { err = terr; OMP("omp flush(err)"); }
	if (terr || bfs_tree[i] < 0 /* both i & j have the same sign */)
	  continue;

	/* Both i and j are in the tree, count as a traversed edge.

	   NOTE: This counts self-edges and repeated edges.  They're
	   part of the input data.
	*/

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
	  terr = -14;
	if (terr) { err = terr; OMP("omp flush(err)"); }
      }

    if (!terr) {
      /* Check that every BFS edge was seen and that there's only one root. */
      OMP("omp for") MTA("mta assert parallel") MTA("mta use 100 streams")
	for (k = 0; k < NV; ++k) {
	  if (level[k] > tmaxdepth) tmaxdepth = level[k];
	  terr = err;
	  if (!terr && k != root) {
	    if (bfs_tree[k] >= 0 && !seen_edge[k])
	      terr = -15;
	    if (bfs_tree[k] == k)
	      terr = -16;
	    if (terr) { err = terr; OMP("omp flush(err)"); }
	  }
	}
    }

    OMP("omp critical") {
      if (tmaxdepth > maxdepth) maxdepth = tmaxdepth;
    }
  }
 done:

  xfree_large (seen_edge);
  if (err) return err;
  return maxdepth;
}
