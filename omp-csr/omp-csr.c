/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010-2011,  Georgia Institute of Technology, USA. */
/* Copyright 2013,  Regents of the University of California, USA. */
/* See COPYING for license. */
#include "../compat.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <assert.h>

#include <alloca.h>

static int64_t int64_fetch_add (int64_t* p, int64_t incr);
static int64_t int64_casval(int64_t* p, int64_t oldval, int64_t newval);
static int int64_cas(int64_t* p, int64_t oldval, int64_t newval);

#include "../graph500.h"
#include "../xalloc.h"
#include "../packed_edge.h"
#include "../timer.h"
#include "../generator.h"
#include "../sorts.h"

#include "bitmap.h"

#define THREAD_BUF_LEN 16384
#define ALPHA 14
#define BETA  24

char IMPLEMENTATION[] = "Reference OpenMP";

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

static int64_t
prefix_sum (const size_t N, int64_t * restrict A, int64_t * restrict buf)
{
  int nt, tid;
  int64_t slice_begin, slice_end, t1, t2, k;

  nt = omp_get_num_threads ();
  tid = omp_get_thread_num ();

  t1 = N / nt;
  t2 = N % nt;
  slice_begin = t1 * tid + (tid < t2? tid : t2);
  slice_end = t1 * (tid+1) + ((tid+1) < t2? (tid+1) : t2);

  buf[tid] = 0;
  for (k = slice_begin; k < slice_end; ++k)
    buf[tid] += A[k];
  OMP("omp barrier");
  OMP("omp single")
    for (k = 1; k < nt; ++k)
      buf[k] += buf[k-1];
  if (tid)
    t1 = buf[tid-1];
  else
    t1 = 0;
  for (k = slice_begin; k < slice_end; ++k) {
    int64_t tmp = A[k];
    A[k] = t1;
    t1 += tmp;
  }
  OMP("omp flush (A)");
  OMP("omp barrier");
  return buf[nt-1];
}

int
create_graph_from_edgelist (struct packed_edge *IJ, int64_t nedge, int64_t nv_in)
{
  int64_t * restrict xoff_half = NULL;
  int64_t * restrict xadj_half = NULL;
  int64_t *buf = NULL;
  size_t sz;
  int64_t accum;
  int errflag = 0;

  int64_t nself_edges = 0;
  int64_t ndup = 0;

  nv = nv_in;
  maxvtx = nv-1;
  xoff = NULL;
  xadj = NULL;

  /* Allocate offset arrays.  */
  sz = (nv+1) * sizeof (*xoff);
  xoff = xmalloc_large_ext (sz);
  xoff_half = xmalloc_large_ext (sz);
  if (!xoff || !xoff_half) { errflag = 1; goto end; }

  OMP("omp parallel") {
    OMP("omp single")
      buf = xmalloc (omp_get_num_threads () * sizeof (*buf));
    if (!buf) { errflag = 1; goto inner_err; }

    /* Set up xoff_half to hold a single copy of the edges. */
    OMP("omp for")
      for (int64_t k = 0; k <= nv; ++k)
	xoff_half[k] = 0;

    OMP("omp for reduction(+:nself_edges)")
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
	  OMP("omp atomic")
	    ++xoff_half[i+1];
	} else
	  ++nself_edges;
      }

    /* Prefix sum to convert counts to offsets. */
    accum = prefix_sum (nv, &xoff_half[1], buf);

    assert (accum + nself_edges == nedge);
    /* All either counted or ignored. */

    /* Allocate endpoint storage.  Do not yet know final number of edges. */
    OMP("omp single")
      xadj_half = xmalloc_large_ext (accum * sizeof (*xadj_half));
    if (!xadj_half) { errflag = 1; goto inner_err; }

    /* Copy endpoints. */
    OMP("omp for")
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
	  where = int64_fetch_add (&xoff_half[i+1], 1);
	  xadj_half[where] = j;
	}
      }

    OMP("omp for")
      for (int64_t i = 0; i <= nv; ++i)
	xoff[i] = 0;

    /* Collapse duplicates and count final numbers. */
    OMP("omp for reduction(+:ndup)")
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
	OMP("omp atomic")
	  xoff[i+1] += kcur - xoff_half[i];

	/* Scatter adjacent vertex j counts. */
	for (int64_t k = xoff_half[i]; k < kcur; ++k) {
	  const int64_t j = xadj_half[k];
	  OMP("omp atomic") ++xoff[j+1];
	}
      }

    assert (xoff[0] == 0);

    /* Another prefix sum to convert counts to offsets. */
    accum = prefix_sum (nv, &xoff[1], buf);

    assert (xoff[0] == 0);

    assert (accum % 2 == 0);
    if (accum/2 + ndup + nself_edges != nedge) {
      fprintf (stderr, "%" PRId64" / 2 + %" PRId64 " + %" PRId64 " (%" PRId64 ") != %" PRId64 "   ( %" PRId64 " )\n",
	       accum, ndup, nself_edges, (accum/2 + ndup + nself_edges), nedge,
	       nedge - (accum/2 + ndup + nself_edges) );
    }
    assert (accum/2 + ndup + nself_edges == nedge);

    /* Allocate final endpoint storage. */
    OMP("omp single")
      xadj = xmalloc_large_ext (accum * sizeof (*xadj));
    if (!xadj) { errflag = 1; goto inner_err; }

    /* Copy to final locations. */
    OMP("omp for")
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
      }

    OMP("omp for")
      for (int64_t i = 0; i < nv; ++i) {
	/* Scatter the other side. */
	const int64_t khalfstart = xoff_half[i];
	const int64_t khalfend = xoff_half[i+1];
	for (int64_t k = khalfstart; k < khalfend; ++k) {
	  const int64_t j = xadj_half[k];
	  if (j < 0) break;
	  const int64_t where = int64_fetch_add (&xoff[j+1], 1);
	  xadj[where] = i;
	}
      }

    assert (accum == xoff[nv]);

  inner_err:
    /* Errors handled after threads join. */;
  }

 end:
  if (buf) free (buf);
  if (xoff_half) xfree_large (xoff_half);
  if (xadj_half) xfree_large (xadj_half);
  if (errflag) {
    if (xoff) xfree_large (xoff);
    if (xadj) xfree_large (xadj);
    xoff = NULL;
    xadj = NULL;
  }
  return errflag;
}

static void
fill_bitmap_from_queue(bitmap_t *bm, int64_t *vlist, int64_t out, int64_t in)
{
  OMP("omp for")
    for (long q_index=out; q_index<in; q_index++)
      bm_set_bit_atomic(bm, vlist[q_index]);
}

static void
fill_queue_from_bitmap(bitmap_t *bm, int64_t *vlist, int64_t *out, int64_t *in,
		       int64_t *local)
{
  OMP("omp single") {
    *out = 0;
    *in = 0;
  }
  OMP("omp barrier");
  int64_t nodes_per_thread = (nv + omp_get_num_threads() - 1) /
			     omp_get_num_threads();
  int64_t i = nodes_per_thread * omp_get_thread_num();
  int local_index = 0;
  if (i < nv) {
    int64_t i_end = i + nodes_per_thread;
    if (i_end >= nv)
      i_end = nv;
    if (bm_get_bit(bm, i))
      local[local_index++] = i;
    i = bm_get_next_bit(bm,i);
    while ((i != -1) && (i < i_end)) {
      local[local_index++] = i;
      i = bm_get_next_bit(bm,i);
      if (local_index == THREAD_BUF_LEN) {
	int my_in = int64_fetch_add(in, THREAD_BUF_LEN);
	for (local_index=0; local_index<THREAD_BUF_LEN; local_index++) {
	  vlist[my_in + local_index] = local[local_index];
	}
	local_index = 0;
      }
    }
  }
  int my_in = int64_fetch_add(in, local_index);
  for (int i=0; i<local_index; i++) {
    vlist[my_in + i] = local[i];
  }
}

static int64_t
bfs_bottom_up_step(int64_t *bfs_tree, bitmap_t *past, bitmap_t *next)
{
  OMP("omp single") {
    bm_swap(past, next);
  }
  OMP("omp barrier");
  bm_reset(next);
  static int64_t awake_count;
  OMP("omp single")
    awake_count = 0;
  OMP("omp barrier");
  OMP("omp for reduction(+ : awake_count)")
    for (int64_t i=0; i<nv; i++) {
      if (bfs_tree[i] == -1) {
	  for (int64_t vo = xoff[i]; vo < xoff[1+i]; vo++) {
	    const int64_t j = xadj[vo];
	  if (bm_get_bit(past, j)) {
	    // printf("%lu\n",i);
	    bfs_tree[i] = j;
	    bm_set_bit_atomic(next, i);
	    awake_count++;
	    break;
	  }
	}
      }
    }
  OMP("omp barrier");
  return awake_count;
}

static void
bfs_top_down_step(int64_t *bfs_tree, int64_t *vlist, int64_t *local, int64_t *k1_p, int64_t *k2_p)
{
  const int64_t oldk2 = *k2_p;
  int64_t kbuf = 0;
  OMP("omp barrier");
  OMP("omp for")
	for (int64_t k = *k1_p; k < oldk2; ++k) {
	  const int64_t v = vlist[k];
	  const int64_t veo = xoff[1+v];
	  int64_t vo;
	  for (vo = xoff[v]; vo < veo; ++vo) {
	    const int64_t j = xadj[vo];
	    if (bfs_tree[j] == -1) {
	      if (int64_cas (&bfs_tree[j], -1, v)) {
		if (kbuf < THREAD_BUF_LEN) {
		  local[kbuf++] = j;
		} else {
		  int64_t voff = int64_fetch_add (k2_p, THREAD_BUF_LEN), vk;
		  assert (voff + THREAD_BUF_LEN <= nv);
		  for (vk = 0; vk < THREAD_BUF_LEN; ++vk)
		    vlist[voff + vk] = local[vk];
		  local[0] = j;
		  kbuf = 1;
		}
	      }
	    }
	  }
	}
  if (kbuf) {
	int64_t voff = int64_fetch_add (k2_p, kbuf), vk;
	assert (voff + kbuf <= nv);
	for (vk = 0; vk < kbuf; ++vk)
	  vlist[voff + vk] = local[vk];
  }
  OMP("omp single")
    *k1_p = oldk2;
  OMP("omp barrier");
  return;
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

  vlist[0] = srcvtx;
  k1 = 0; k2 = 1;
  bfs_tree[srcvtx] = srcvtx;

  bitmap_t past, next;
  bm_init(&past, nv);
  bm_init(&next, nv);

  int64_t down_cutoff = nv / BETA;
  int64_t scout_count = xoff[1+srcvtx] - xoff[srcvtx];

  OMP("omp parallel shared(k1, k2, scout_count)") {
    int64_t k;
    int64_t nbuf[THREAD_BUF_LEN];
    int64_t awake_count = 1;
    int64_t edges_to_check = xoff[nv];

    OMP("omp for")
      for (k = 0; k < srcvtx; ++k)
	bfs_tree[k] = -1;
    OMP("omp for")
      for (k = srcvtx+1; k < nv; ++k)
	bfs_tree[k] = -1;

    while (awake_count != 0) {
      // Top-down
      if (scout_count < ((edges_to_check - scout_count)/ALPHA)) {
	bfs_top_down_step(bfs_tree, vlist, nbuf, &k1, &k2);
	edges_to_check -= scout_count;
	awake_count = k2-k1;
      // Bottom-up
      } else {
	fill_bitmap_from_queue(&next, vlist, k1, k2);
	do {
	  awake_count = bfs_bottom_up_step(bfs_tree, &past, &next);
	} while ((awake_count > down_cutoff));
	fill_queue_from_bitmap(&next, vlist, &k1, &k2, nbuf);
	OMP("omp barrier");
      }
      // Count the number of edges in the frontier
      OMP("omp single")
	scout_count = 0;
      OMP("omp for reduction(+ : scout_count)")
      for (int64_t i=k1; i<k2; i++) {
	int64_t v = vlist[i];
	scout_count += xoff[1+v] - xoff[v];
      }
    }
  }

  bm_free(&past);
  bm_free(&next);
  xfree_large (vlist);

  return err;
}

void
destroy_graph (void)
{
  free_graph ();
}

#if defined(_OPENMP)
#if defined(__GNUC__)||defined(__INTEL_COMPILER)
int64_t
int64_fetch_add (int64_t* p, int64_t incr)
{
  return __sync_fetch_and_add (p, incr);
}
int64_t
int64_casval(int64_t* p, int64_t oldval, int64_t newval)
{
  return __sync_val_compare_and_swap (p, oldval, newval);
}
int
int64_cas(int64_t* p, int64_t oldval, int64_t newval)
{
  return __sync_bool_compare_and_swap (p, oldval, newval);
}
#else
/* XXX: These are not correct, but suffice for the above uses. */
int64_t
int64_fetch_add (int64_t* p, int64_t incr)
{
  int64_t t;
  OMP("omp critical") {
    t = *p;
    *p += incr;
  }
  OMP("omp flush (p)");
  return t;
}
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
int
int64_cas(int64_t* p, int64_t oldval, int64_t newval)
{
  int out = 0;
  OMP("omp critical (CAS)") {
    int64_t v = *p;
    if (v == oldval) {
      *p = newval;
      out = 1;
    }
  }
  OMP("omp flush (p)");
  return out;
}
#endif
#else
int64_t
int64_fetch_add (int64_t* p, int64_t incr)
{
  int64_t t = *p;
  *p += incr;
  return t;
}
int64_t
int64_casval(int64_t* p, int64_t oldval, int64_t newval)
{
  int64_t v = *p;
  if (v == oldval)
    *p = newval;
  return v;
}
int
int64_cas(int64_t* p, int64_t oldval, int64_t newval)
{
  int64_t v = *p;
  int out = 0;
  if (v == oldval) {
    *p = newval;
    out = 1;
  }
  return out;
}
#endif
