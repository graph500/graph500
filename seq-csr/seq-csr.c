/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010-2011,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#define _FILE_OFFSET_BITS 64
#define _THREAD_SAFE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include "../compat.h"
#include "../graph500.h"
#include "../xalloc.h"
#include "../packed_edge.h"
#include "../generator.h"
#include "../sorts.h"

char IMPLEMENTATION[] = "Reference sequential";

static int64_t maxvtx, nv;
static int64_t * restrict xoff; /* Length 2*nv+2 */
static int64_t * restrict xadj;

static void
free_graph (void)
{
  if (xadj) xfree_large (xadj);
  if (xoff) xfree_large (xoff);
}

static inline void
canonical_order_edge (int64_t * restrict ip, int64_t * restrict jp)
{
  int64_t i = *ip, j = *jp;
#if !defined(ORDERED_TRIL)
  if (((i ^ j) & 0x01) == (i < j)) {
    *ip = i;
    *jp = j;
  } else {
    *ip = j;
    *jp = i;
  }
#else
  if (i < j) {
    *ip = j;
    *jp = i;
  }
#endif
}

int 
create_graph_from_edgelist (struct packed_edge *IJ, int64_t nedge, int64_t nv_in)
{
  int64_t * restrict xoff_half = NULL;
  int64_t * restrict xadj_half = NULL;
  size_t sz;
  int64_t accum;

  nv = nv_in;
  maxvtx = nv-1;
  xoff = NULL;
  xadj = NULL;

  /* Allocate offset arrays.  */
  sz = (nv+1) * sizeof (*xoff);
  xoff = xmalloc_large_ext (sz);
  xoff_half = xmalloc_large_ext (sz);
  if (!xoff || !xoff_half) goto err;

  /* Set up xoff_half to hold a single copy of the edges. */
  for (int64_t k = 0; k <= nv; ++k)
    xoff_half[k] = 0;

  int64_t nself_edges = 0;

  for (int64_t k = 0; k < nedge; ++k) {
    int64_t i, j;
#if !defined(STORED_EDGELIST)
    uint8_t w;
    make_edge (k, &i, &j, &w);
#else
    i = get_v0_from_edge(&IJ[k]);
    j = get_v1_from_edge(&IJ[k]);
#endif
    assert (i >= 0);
    assert (j >= 0);
    if (i != j) { /* Skip self-edges. */
      canonical_order_edge (&i, &j);
      ++xoff_half[i+1];
    } else
      ++nself_edges;
  }

  /* Prefix sum to convert counts to offsets. */
  accum = 0;
  for (int64_t k = 1; k <= nv; ++k) {
    int64_t tmp = xoff_half[k];
    xoff_half[k] = accum;
    accum += tmp;
  }

  assert (accum + nself_edges == nedge);
  /* All either counted or ignored. */

  /* Allocate endpoint storage.  Do not yet know final number of edges. */
  xadj_half = xmalloc_large_ext (accum * sizeof (*xadj_half));
  if (!xadj_half) goto err;

  /* Copy endpoints. */
  for (int64_t k = 0; k < nedge; ++k) {
    int64_t i, j;
#if !defined(STORED_EDGELIST)
    uint8_t w;
    make_edge (k, &i, &j, &w);
#else
    i = get_v0_from_edge(&IJ[k]);
    j = get_v1_from_edge(&IJ[k]);
#endif
    assert (i >= 0);
    assert (j >= 0);
    if (i != j) { /* Skip self-edges. */
      int64_t where;
      canonical_order_edge (&i, &j);
      where = xoff_half[i+1]++;
      xadj_half[where] = j;
    }
  }

  int64_t ndup = 0;

  for (int64_t i = 0; i <= nv; ++i)
    xoff[i] = 0;

  /* Collapse duplicates and count final numbers. */
  for (int64_t i = 0; i < nv; ++i) {
    int64_t kcur = xoff_half[i];
    const int64_t kend = xoff_half[i+1];

    if (kcur == kend) continue; /* Empty. */

    introsort_i64 (&xadj_half[kcur], kend-kcur);

    for (int64_t k = kcur+1; k < kend; ++k)
      if (xadj_half[k] != xadj_half[kcur]) {
	xadj_half[++kcur] = xadj_half[k];
      }
    ++kcur;
    if (kcur != kend)
      xadj_half[kcur] = -1;
    ndup += kend - kcur;

    /* Current vertex i count. */
    xoff[i+1] += kcur - xoff_half[i];

    /* Scatter adjacent vertex j counts. */
    for (int64_t k = xoff_half[i]; k < kcur; ++k) {
      const int64_t j = xadj_half[k];
      ++xoff[j+1];
    }
  }

  assert (xoff[0] == 0);

  /* Another prefix sum to convert counts to offsets. */
  accum = 0;
  for (int64_t k = 1; k <= nv; ++k) {
    int64_t tmp = xoff[k];
    xoff[k] = accum;
    accum += tmp;
  }

  assert (xoff[0] == 0);

  assert (accum % 2 == 0);
  if (accum/2 + ndup + nself_edges != nedge) {
    fprintf (stderr, "%" PRId64" / 2 + %" PRId64 " + %" PRId64 " (%" PRId64 ") != %" PRId64 "   ( %" PRId64 " )\n",
	     accum, ndup, nself_edges, (accum/2 + ndup + nself_edges), nedge,
	     nedge - (accum/2 + ndup + nself_edges) );
  }
  assert (accum/2 + ndup + nself_edges == nedge);

  /* Allocate final endpoint storage. */
  xadj = xmalloc_large_ext (accum * sizeof (*xadj));
  if (!xadj) goto err;

  /* Copy to final locations. */
  for (int64_t i = 0; i < nv; ++i) {
    const int64_t khalfstart = xoff_half[i];
    const int64_t khalfend = xoff_half[i+1];
    int64_t khalflen = khalfend - khalfstart;
    const int64_t kstart = xoff[i+1];
    int64_t k;

    if (!khalflen) continue; /* Empty. */

    /* Copy the contiguous portion. */
    for (k = 0; k < khalflen; ++k) {
      const int64_t j = xadj_half[khalfstart + k];
      if (j < 0) break;
      xadj[kstart + k] = j;
    }
    xoff[i+1] += k;

    /* Scatter the other side. */
    for (int64_t k2 = 0; k2 < k; ++k2) {
      const int64_t j = xadj_half[khalfstart + k2];
      const int64_t where = xoff[j+1]++;
      xadj[where] = i;
    }
  }

  assert (accum == xoff[nv]);

  xfree_large (xoff_half);
  xfree_large (xadj_half);
  return 0;

 err:
  if (xoff_half) xfree_large (xoff_half);
  if (xadj_half) xfree_large (xadj_half);
  if (xoff) xfree_large (xoff);
  if (xadj) xfree_large (xadj);
  return -1;
}

int
make_bfs_tree (int64_t *bfs_tree_out, int64_t *max_vtx_out,
	       int64_t srcvtx)
{
  int64_t * restrict bfs_tree = bfs_tree_out;
  int err = 0;

  int64_t * restrict vlist = NULL;
  int64_t k1, k2;

  *max_vtx_out = maxvtx;

  vlist = xmalloc_large (nv * sizeof (*vlist));
  if (!vlist) return -1;

  for (k1 = 0; k1 < nv; ++k1)
    bfs_tree[k1] = -1;

  vlist[0] = srcvtx;
  bfs_tree[srcvtx] = srcvtx;
  k1 = 0; k2 = 1;
  while (k1 != k2) {
    const int64_t oldk2 = k2;
    int64_t k;
    for (k = k1; k < oldk2; ++k) {
      const int64_t v = vlist[k];
      const int64_t veo = xoff[v+1];
      int64_t vo;
      for (vo = xoff[v]; vo < veo; ++vo) {
	const int64_t j = xadj[vo];
	if (bfs_tree[j] == -1) {
	  bfs_tree[j] = v;
	  vlist[k2++] = j;
	}
      }
    }
    k1 = oldk2;
  }

  xfree_large (vlist);

  return err;
}

void
destroy_graph (void)
{
  free_graph ();
}
