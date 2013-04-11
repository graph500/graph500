/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010-2011,  Georgia Institute of Technology, USA. */
/* Copyright 2013,  Regents of the University of California, USA. */
/* See COPYING for license. */
#include "../compat.h"
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <assert.h>

static int64_t int64_fetch_add (int64_t* p, int64_t incr);
static uint64_t uint64_fetch_or (uint64_t* p, uint64_t v);
static int int64_cas(int64_t* p, int64_t oldval, int64_t newval);

#include "../graph500.h"
#include "../xalloc.h"
#include "../packed_edge.h"
#include "../timer.h"
#include "../generator.h"
#include "../sorts.h"

#define THREAD_BUF_LEN 256
#define ALPHA 14
#define BETA  24

#define WORD_OFFSET(n) (n/64)
#define BIT_OFFSET(n) (n & 0x3f)

char IMPLEMENTATION[] = "Reference OpenMP";

static int64_t nv;
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

static inline void
flush_buf (const int64_t * restrict buf, const int blen,
	  int64_t * restrict glist, int64_t * restrict glen)
{
  if (blen) {
    const int64_t base = int64_fetch_add (glen, blen);
    for (int k = 0; k < blen; ++k)
      glist[base + k] = buf[k];
  }
}

static inline void
push_buf (int64_t * restrict buf, int * restrict blen,
	  int64_t k,
	  int64_t * restrict glist, int64_t * restrict glen)
{
  int64_t where;
  if (*blen == THREAD_BUF_LEN) {
    flush_buf (buf, THREAD_BUF_LEN, glist, glen);
    *blen = 0;
  }
  where = *blen;
  ++*blen;
  buf[where] = k;
}

#define BIT(i) (((uint64_t)1)<<BIT_OFFSET((i)))
#define ATOMIC_SET_BIT(bm, i) (uint64_fetch_or (&((bm)[WORD_OFFSET((i))]), BIT(i)))
#define FETCH_SET_BIT(bm, i) (ATOMIC_SET_BIT(bm, i) & BIT(i))
#define SET_BIT(bm, i) ((bm)[WORD_OFFSET((i))] |= BIT(i))
#define GET_BIT(bm, i) ((bm)[WORD_OFFSET((i))] & BIT(i))

static void
fill_bitmap_from_queue (uint64_t * restrict bm, size_t bmlen,
			int64_t * restrict vlist,
			int64_t out, int64_t in)
{
  OMP("omp for")
    for (int64_t k = 0; k < bmlen; ++k)
      bm[k] = 0;
  OMP("omp for")
    for (int64_t q_index = out; q_index < in; q_index++)
      ATOMIC_SET_BIT (bm, vlist[q_index]);
}

static void
fill_queue_from_bitmap (uint64_t * restrict bm,
			int64_t * restrict vlist,
			int64_t * restrict head, int64_t * restrict tail,
			int64_t * restrict local)
{
  int len = 0;

  OMP("omp single") {
    *head = 0;
    *tail = 0;
  }

  OMP("omp for nowait")
    for (int64_t k = 0; k < nv; k += 64) {
      const uint64_t m = bm[k/64];
      if (m)
	for (int kk = 0; kk < 64; ++kk)
	  if (m & (((uint64_t)1) << kk))
	    push_buf (local, &len, k+kk, vlist, tail);
    }
  flush_buf (local, len, vlist, tail);
  OMP("omp barrier");
}

static void
bfs_bottom_up_step (int64_t * restrict bfs_tree,
		    const uint64_t * restrict past, uint64_t * restrict next,
		    const size_t bmlen,
		    int64_t * awake_count)
{
  int64_t t_awake_count = 0;

  OMP("omp for")
    for (int64_t k = 0; k < bmlen; ++k)
      next[k] = 0;

  OMP("omp for nowait")
    for (int64_t i = 0; i < nv; i++) {
      if (bfs_tree[i] == -1) {
	for (int64_t vo = xoff[i]; vo < xoff[1+i]; vo++) {
	  const int64_t j = xadj[vo];
	  if (GET_BIT (past, j)) {
	    ATOMIC_SET_BIT (next, i);
	    bfs_tree[i] = j;
	    t_awake_count++;
	    break;
	  }
	}
      }
    }

  OMP("omp atomic") *awake_count += t_awake_count;
  OMP("omp barrier");
}

static void
bfs_top_down_step (int64_t * restrict bfs_tree,
		   int64_t * restrict vlist,
		   int64_t * head_in, int64_t * tail_in,
		   int64_t * restrict local)
{
  const int64_t tail = *tail_in;
  const int64_t head = *head_in;
  int len = 0;

  OMP("omp barrier");
  OMP("omp for nowait")
    for (int64_t k = head; k < tail; ++k) {
      const int64_t v = vlist[k];
      const int64_t veo = xoff[1+v];
      int64_t vo;
      for (vo = xoff[v]; vo < veo; ++vo) {
	const int64_t j = xadj[vo];
	if (bfs_tree[j] == -1)
	  if (int64_cas (&bfs_tree[j], -1, v))
	    push_buf (local, &len, j, vlist, tail_in);
      }
    }
  flush_buf (local, len, vlist, tail_in);

  OMP("omp single")
    *head_in = tail;

  return;
}

int
make_bfs_tree (int64_t *bfs_tree_out, int64_t srcvtx)
{
  int64_t * restrict bfs_tree = bfs_tree_out;
  const size_t bmlen = (nv + 63) & ~63;
  uint64_t * restrict past = NULL;
  uint64_t * restrict next = NULL;
  int err = 0;

  int64_t * restrict vlist = NULL;
  int64_t k1, k2;

  const int64_t down_cutoff = nv / BETA;
  int64_t scout_count = xoff[1+srcvtx] - xoff[srcvtx];
  int64_t awake_count = 1;

  /* Size sanity checks. */
#if 8 != CHAR_BIT
#error "Hard-coded to support eight-bit bytes."
#endif
  assert (8 == sizeof (uint64_t));

  vlist = xmalloc_large (nv * sizeof (*vlist));
  if (!vlist) return -1;

  past = xmalloc_large (bmlen * sizeof (*past));
  next = xmalloc_large (bmlen * sizeof (*next));
  if (!past || !next) { err = -1; goto done; }

  vlist[0] = srcvtx;
  k1 = 0; k2 = 1;

  OMP("omp parallel") {
#if defined(MALLOC_THREAD_BUF)
    int64_t * restrict local = xmalloc (THREAD_BUF_LEN * sizeof (*local));
#else
    int64_t local[THREAD_BUF_LEN];
#endif
    int64_t edges_to_check = xoff[nv];

    OMP("omp for")
      for (int64_t k = 0; k < nv; ++k)
	bfs_tree[k] = (k == srcvtx? srcvtx : -1);

    while (k1 != k2) {
      OMP("omp barrier");
      if (scout_count < ((edges_to_check - scout_count)/ALPHA)) {
	// Top-down
	bfs_top_down_step (bfs_tree, vlist, &k1, &k2, local);
	edges_to_check -= scout_count;
      } else {
	// Bottom-up
	fill_bitmap_from_queue (next, bmlen, vlist, k1, k2);
	do {
	  OMP("omp barrier");
	  /* Swap past & next */
	  OMP("omp single") {
	    uint64_t * restrict t = past;
	    past = next;
	    next = t;
	    awake_count = 0;
	  }
	  bfs_bottom_up_step (bfs_tree, past, next, bmlen, &awake_count);
	} while (awake_count > down_cutoff);
	fill_queue_from_bitmap (next, vlist, &k1, &k2, local);
	assert (awake_count == (k2 - k1));
      }
      // Count the number of edges in the frontier
      OMP("omp single")
	scout_count = 0;
      OMP("omp for reduction(+: scout_count)")
	for (int64_t i=k1; i<k2; i++) {
	  int64_t v = vlist[i];
	  scout_count += xoff[1+v] - xoff[v];
	}
    }

#if defined(MALLOC_THREAD_BUF)
    free (local);
#endif
  }

 done:
  if (past) xfree_large (past);
  if (next) xfree_large (next);
  if (vlist) xfree_large (vlist);

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
uint64_t
uint64_fetch_or (uint64_t* p, uint64_t v)
{
  return __sync_fetch_and_or (p, v);
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
uint64_t
uint64_fetch_or (uint64_t* p, uint64_t v)
{
  uint64_t t;
  OMP("omp critical") {
    t = *p;
    *p |= v;
  }
  OMP("omp flush (p)");
  return t;
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
int64_fetch_or (int64_t* p, int64_t incr)
{
  int64_t t = *p;
  *p |= incr;
  return t;
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
