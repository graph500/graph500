/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
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
#include "../generator/graph_generator.h"

#define MINVECT_SIZE 2

static int64_t maxvtx, nv, sz;
static int64_t * restrict xoff; /* Length 2*nv+2 */
static int64_t * restrict xadjstore; /* Length MINVECT_SIZE + (xoff[nv] == nedge) */
static int64_t * restrict xadj;

static void
find_nv (const struct packed_edge * restrict IJ, const int64_t nedge)
{
  int64_t k, tmaxvtx = -1;

  maxvtx = -1;
  for (k = 0; k < nedge; ++k) {
    if (get_v0_from_edge(&IJ[k]) > tmaxvtx)
      tmaxvtx = get_v0_from_edge(&IJ[k]);
    if (get_v1_from_edge(&IJ[k]) > tmaxvtx)
      tmaxvtx = get_v1_from_edge(&IJ[k]);
  }
  k = readfe (&maxvtx);
  if (tmaxvtx > k)
    k = tmaxvtx;
  writeef (&maxvtx, k);
  nv = 1+maxvtx;
}

static int
alloc_graph (int64_t nedge)
{
  sz = (2*nv+2) * sizeof (*xoff);
  xoff = xmalloc_large (sz);
  if (!xoff) return -1;
  return 0;
}

static void
free_graph (void)
{
  xfree_large (xadjstore);
  xfree_large (xoff);
}

#define XOFF(k) (xoff[2*(k)])
#define XENDOFF(k) (xoff[1+2*(k)])

static int
setup_deg_off (const struct packed_edge * restrict IJ, int64_t nedge)
{
  int64_t k, accum;
  for (k = 0; k < 2*nv+2; ++k)
    xoff[k] = 0;
  MTA("mta assert nodep") MTA("mta use 100 streams")
  for (k = 0; k < nedge; ++k) {
    int64_t i = get_v0_from_edge(&IJ[k]);
    int64_t j = get_v1_from_edge(&IJ[k]);
    if (i != j) { /* Skip self-edges. */
      int_fetch_add (&XOFF(i), 1);
      int_fetch_add (&XOFF(j), 1);
    }
  }
  accum = 0;
  MTA("mta use 100 streams")
  for (k = 0; k < nv; ++k) {
    int64_t tmp = XOFF(k);
    if (tmp < MINVECT_SIZE) tmp = MINVECT_SIZE;
    XOFF(k) = accum;
    accum += tmp;
  }
  XOFF(nv) = accum;
  MTA("mta use 100 streams")
  for (k = 0; k < nv; ++k)
    XENDOFF(k) = XOFF(k);
  if (!(xadjstore = malloc ((accum + MINVECT_SIZE) * sizeof (*xadjstore))))
    return -1;
  xadj = &xadjstore[MINVECT_SIZE]; /* Cheat and permit xadj[-1] to work. */
  for (k = 0; k < accum + MINVECT_SIZE; ++k)
    xadjstore[k] = -1;
  return 0;
}

static void
scatter_edge (const int64_t i, const int64_t j)
{
  int64_t where;
  where = int_fetch_add (&XENDOFF(i), 1);
  xadj[where] = j;
}

static int
i64cmp (const void *a, const void *b)
{
  const int64_t ia = *(const int64_t*)a;
  const int64_t ib = *(const int64_t*)b;
  if (ia < ib) return -1;
  if (ia > ib) return 1;
  return 0;
}

MTA("mta no inline")
static void
pack_edges (void)
{
  int64_t v;

  MTA("mta assert parallel") MTA("mta use 100 streams")
    for (v = 0; v < nv; ++v) {
      int64_t kcur, k;
      if (XOFF(v)+1 < XENDOFF(v)) {
	qsort (&xadj[XOFF(v)], XENDOFF(v)-XOFF(v), sizeof(*xadj), i64cmp);
	kcur = XOFF(v);
	MTA("mta loop serial")
        for (k = XOFF(v)+1; k < XENDOFF(v); ++k)
	  if (xadj[k] != xadj[kcur])
	    xadj[1 + int_fetch_add (&kcur, 1)] = xadj[k];
	++kcur;
	for (k = kcur; k < XENDOFF(v); ++k)
	  xadj[k] = -1;
	XENDOFF(v) = kcur;
      }
    }
}

static void
gather_edges (const struct packed_edge * restrict IJ, int64_t nedge)
{
  int64_t k;

  MTA("mta assert nodep")
  MTA("mta use 100 streams")
  for (k = 0; k < nedge; ++k) {
    int64_t i = get_v0_from_edge(&IJ[k]);
    int64_t j = get_v1_from_edge(&IJ[k]);
    if (i != j) {
      scatter_edge (i, j);
      scatter_edge (j, i);
    }
  }

  pack_edges ();
}

int 
create_graph_from_edgelist (struct packed_edge *IJ, int64_t nedge)
{
  find_nv (IJ, nedge);
  if (alloc_graph (nedge)) return -1;
  if (setup_deg_off (IJ, nedge)) {
    xfree_large (xoff);
    return -1;
  }
  gather_edges (IJ, nedge);
  return 0;
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
    MTA("mta assert nodep") MTA("mta use 100 streams") MTA("mta interleave schedule")
    for (k = k1; k < oldk2; ++k) {
      const int64_t v = vlist[k];
      const int64_t veo = XENDOFF(v);
      int64_t vo;
      MTA("mta loop serial")
      for (vo = XOFF(v); vo < veo; ++vo) {
	const int64_t j = xadj[vo];
	if (bfs_tree[j] == -1) {
	  int64_t t = readfe (&bfs_tree[j]);
	  if (t == -1) {
	    t = v;
	    vlist[int_fetch_add (&k2, 1)] = j;
	  }
	  writeef (&bfs_tree[j], t);
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
