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
#include "prng.h"
#include "generator.h"
#include "verify.h"

static int64_t int64_casval(int64_t* p, int64_t oldval, int64_t newval);

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

static size_t
bsearch_i64 (const int64_t * restrict m, const size_t sz, const int64_t key)
{
  size_t start = 0;
  size_t end = sz;

  if (m[start] == key) return start;
  if (m[end-1] == key) return end-1;

  while (end - start > 8) {
    size_t midpt = (start + end)/2;
    if (m[midpt] == key) return midpt;
    if (m[midpt] > key) end = midpt;
    else start = midpt;
  }
  for (; start < end; ++start)
    if (m[start] == key) return start;
  return (size_t)-1;
}

static void
gather_sample (const int64_t * restrict tree,
	       const int64_t * restrict depth,
	       int64_t * restrict tree_edge,
	       const int64_t * restrict sampled_vertex,
	       int64_t * restrict sampled_edge,
	       const struct packed_edge * restrict IJ)
{
  const size_t NSAMPV = (2*SCALE > NV? NV : 2*SCALE);

  OMP("omp parallel") {
    OMP("omp for nowait")
      for (size_t k = 0; k < NSAMPV*EF; ++k) {
	sampled_edge[k] = -1;
	sampled_edge[NSAMPV*EF + k] = 0;
      }

    OMP("omp for")
      for (int64_t k = 0; k < NV; ++k)
	tree_edge[k] = 0;

    OMP("omp for")
      for (int64_t k = 0; k < NE; ++k) {
	int64_t i, j;
	uint8_t w = 1;

#if !defined(STORED_EDGELIST)
	make_edge (k, &i, &j, &w);
#else
	i = get_v0_from_edge(&IJ[k]);
	j = get_v1_from_edge(&IJ[k]);
#endif

	if (i == j) continue; /* Ignore self edges. */

	/* Orient away from the root. */
	if (depth[i] > depth[j]) {
	  int64_t t = i;
	  i = j;
	  j = t;
	}

	if (tree[j] == i) {
	  /* Tree edge. */
	  OMP("omp atomic")
	    tree_edge[j] += w;
	} else {
	  /* For each of i and j, check if a sampled edge. */
	  size_t i_off = bsearch_i64 (sampled_vertex, NSAMPV,
				      i);
	  size_t j_off = bsearch_i64 (sampled_vertex, NSAMPV,
				      j);

	  if (i_off != (size_t)-1) {
	    size_t base = i_off * EF;
	    for (int64_t kk = 0; kk < EF; ++kk) {
	      int64_t tmp = sampled_edge[base + kk];
	      while (tmp < 0)
		tmp = int64_casval (&sampled_edge[base + kk], tmp, j);
	      if (sampled_edge[kk + i_off*EF] == j) {
		OMP("omp atomic")
		  sampled_edge[kk + i_off*EF + NSAMPV*EF] += w;
		break;
	      }
	    }
	  }

	  if (j_off != (size_t)-1) {
	    size_t base = j_off * EF;
	    for (int64_t kk = 0; kk < EF; ++kk) {
	      int64_t tmp = sampled_edge[base + kk];
	      while (tmp < 0)
		tmp = int64_casval (&sampled_edge[base + kk], tmp, i);
	      if (sampled_edge[base + kk] == i) {
		OMP("omp atomic")
		  sampled_edge[base + kk + NSAMPV*EF] += w;
		break;
	      }
	    }
	  }

	}
      }
  }
}

int64_t
verify_bfs_tree (const int64_t *bfs_tree_in, const int64_t *level_in,
		 const int64_t root, const double kerneltime,
		 const struct packed_edge *IJ_in, int64_t nedge)
{
  const size_t NSAMPV = (2*SCALE > NV? NV : 2*SCALE);
  const int is_bfs = 1;
  const int64_t * restrict bfs_tree = bfs_tree_in;
  const int64_t * restrict level = level_in;
  const struct packed_edge * restrict IJ = IJ_in;

  int err = 0;
  int64_t maxdepth = -1;

  int64_t * restrict level_tmp = NULL;
  int64_t * restrict tree_edge = NULL; /* NV x 2 */
  int64_t * restrict sampled_vertex = NULL; /* NSAMPV */
  int64_t * restrict sampled_edge = NULL; /* NSAMPV x EF x 2*/


  if (root >= NV || bfs_tree[root] != root)
    return -999;

  err = 0;

  OMP("omp parallel") {
    int terr = 0;
    OMP("omp for")
      for (int64_t v = 0; v < NV; ++v) {
	if (bfs_tree[v] < 0 || bfs_tree[v] >= NV)
	  terr = -1000;
      }
    if (terr && !err)
      OMP("omp critical")
	if (!err) err = terr;
  }
  if (err) return err;

  if (!level_in) {
    level_tmp = xmalloc_large (NV * sizeof (*level_tmp));
    err = compute_levels (level_tmp, NV, bfs_tree, root);
    level = level_tmp; /* Type launder to const. */
  }

  if (err) goto done;

  tree_edge = xmalloc_large (2 * NV * sizeof (*tree_edge));
  sampled_vertex = xmalloc (NSAMPV * sizeof (*sampled_vertex));
  sampled_edge = xmalloc (2 * NSAMPV * EF * sizeof (*sampled_edge));

  sample_roots (sampled_vertex, NSAMPV, ceil (NE / kerneltime));
  gather_sample (bfs_tree, level, tree_edge, sampled_vertex, sampled_edge, IJ);

  OMP("omp parallel") {
    int terr = 0;
    int64_t tmaxdepth = 0;
    OMP("omp for")
      for (int64_t k = 0; k < NV; ++k) {
	if (!terr) {
	  /* Check that the tree is rooted and without cycles. */
	  int64_t nhop = 0;
	  int64_t v = k;
	  while (nhop < NV && bfs_tree[v] != v) {
	    v = bfs_tree[v];
	    ++nhop;
	  }
	  if (nhop >= NV) { terr = -1; /* Cycle. */ goto vdone; }
	  if (v != root) { terr = -16; /* More than one root. */ goto vdone; }

	  if (level[k] > tmaxdepth) tmaxdepth = level[k];

	  if (k != root) {
	    /* Check that the tree edge exists. */
	    if (tree_edge[k] == 0)
	      terr = -15;
	    else {
	      /* Check complimentary slackness on the tree edge. */
	      int64_t pd_gap = (is_bfs? 1 : tree_edge[k]) +
		level[bfs_tree[k]] - level[k];
	      if (pd_gap) { terr = -32; goto vdone; }
	    }
	  }
	vdone:
	  ;
	}
      }
    if (terr && !err)
      OMP("omp critical")
	if (!err) err = terr;
    if (tmaxdepth > maxdepth)
      OMP("omp critical")
	if (tmaxdepth > maxdepth) maxdepth = tmaxdepth;
  }

  if (err) goto done; /* Tree itself failed. */

  /* Check sampled edges. */
  OMP("omp parallel") {
    int terr = 0;
    OMP("omp for")
      for (int64_t k = 0; k < NSAMPV; ++k) {
	int64_t i = sampled_vertex[k];
	for (int64_t kk = 0; kk < EF; ++kk) {
	  int64_t j = sampled_edge[kk + k*EF];
	  if (j >= 0) {
	    int64_t w = sampled_edge[kk + k*EF + NSAMPV*EF];
	    int dir = (1 - 2*(level[i] > level[j]));
	    int64_t pd_gap = (is_bfs? 1 : w) + dir * (level[i] - level[j]);
	    if (is_bfs && (pd_gap < -1 || pd_gap > 1))
	      terr = -33; /* Edges cross more than one level. */
	    else if (pd_gap < 0)
	      terr = -34; /* Constraints violated. */
	    if (terr)
	      fprintf (stderr, "rt %d: (%d, %d; %d) [%d]  level[i] %d  level[j] %d  pd_gap %d     %d\n",
		       (int)root, (int)i, (int)j, (int)w, dir, (int)level[i], (int)level[j],
		       (int)pd_gap, terr);
	    assert (!terr);
	  }
	}
      }
    if (terr && !err)
      OMP("omp critical")
	if (!err) err = terr;
  }

 done:

  if (sampled_edge) free (sampled_edge);
  if (sampled_vertex) free (sampled_vertex);
  if (tree_edge) xfree_large (tree_edge);
  if (level_tmp) xfree_large (level_tmp);
  if (err) return err;
  return maxdepth;
}

#if defined(_OPENMP)
#if defined(__GNUC__)||defined(__INTEL_COMPILER)
int64_t
int64_casval(int64_t* p, int64_t oldval, int64_t newval)
{
  return __sync_val_compare_and_swap (p, oldval, newval);
}
#else
/* XXX: These are not correct, but suffice for the above uses. */
int64_t
int64_casval(int64_t* p, int64_t oldval, int64_t newval)
{
  int64_t v;
  OMP("omp critical (CAS)") {
    v = *p;
    if (v == oldval)
      *p = newval;
  }
  OMP("omp flush (p)");
  return v;
}
#endif
#else
int64_t
int64_casval(int64_t* p, int64_t oldval, int64_t newval)
{
  int64_t v = *p;
  if (v == oldval)
    *p = newval;
  return v;
}
#endif
