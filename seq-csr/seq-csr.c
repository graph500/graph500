/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010-2011,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#define _FILE_OFFSET_BITS 64
#define _THREAD_SAFE
#include <stdlib.h>
#include <limits.h>
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
#define ALPHA 14
#define BETA  24

#define WORD_OFFSET(n) (n/64)
#define BIT_OFFSET(n) (n & 0x3f)

static int64_t nv;
static int64_t * restrict xoff; /* Length 2*nv+2 */
static int64_t * restrict xadj;
static int16_t * restrict xval;

static void
free_graph (void)
{
  if (xval) xfree_large (xval);
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
  int16_t * restrict xval_half = NULL;
  size_t sz;
  int64_t accum;

  nv = nv_in;
  xoff = NULL;
  xadj = NULL;
  xval = NULL;

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
    make_edge_endpoints (k, &i, &j);
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
  xval_half = xmalloc_large_ext (accum * sizeof (*xval_half));
  if (!xadj_half || !xval_half) goto err;

  /* Copy endpoints. */
  for (int64_t k = 0; k < nedge; ++k) {
    int64_t i, j;
    uint8_t w;
#if !defined(STORED_EDGELIST)
    make_edge (k, &i, &j, &w);
#else
    i = get_v0_from_edge(&IJ[k]);
    j = get_v1_from_edge(&IJ[k]);
    w = get_w_from_edge(&IJ[k]);
#endif
    assert (i >= 0);
    assert (j >= 0);
    if (i != j) { /* Skip self-edges. */
      int64_t where;
      canonical_order_edge (&i, &j);
      where = xoff_half[i+1]++;
      xadj_half[where] = j;
      xval_half[where] = w;
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

    introsort_both_i64 (&xadj_half[kcur], &xval_half[kcur], kend-kcur);

    for (int64_t k = kcur+1; k < kend; ++k) {
      if (xadj_half[k] != xadj_half[kcur]) {
	xadj_half[++kcur] = xadj_half[k];
	xval_half[kcur] = xval_half[k];
      } else
	xval_half[kcur] += xval_half[k];
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
  assert (accum/2 + ndup + nself_edges == nedge);

  /* Allocate final endpoint storage. */
  xadj = xmalloc_large_ext (accum * sizeof (*xadj));
  xval = xmalloc_large_ext (accum * sizeof (*xval));
  if (!xadj || !xval) goto err;

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
      xval[kstart + k] = xval_half[khalfstart + k];
    }
    xoff[i+1] += k;

    /* Scatter the other side. */
    for (int64_t k2 = 0; k2 < k; ++k2) {
      const int64_t j = xadj_half[khalfstart + k2];
      const int64_t where = xoff[j+1]++;
      xadj[where] = i;
      xval[where] = xval_half[khalfstart + k2];
    }
  }

  assert (accum == xoff[nv]);

  xfree_large (xval_half);
  xfree_large (xadj_half);
  xfree_large (xoff_half);
  return 0;

 err:
  if (xoff_half) xfree_large (xoff_half);
  if (xadj_half) xfree_large (xadj_half);
  if (xval_half) xfree_large (xval_half);
  if (xoff) xfree_large (xoff);
  if (xadj) xfree_large (xadj);
  if (xval) xfree_large (xval);
  xoff = NULL;
  xadj = NULL;
  xval = NULL;
  return -1;
}

#define SET_BIT(bm, i) (bm)[WORD_OFFSET((i))] |= (((uint64_t) 1)<<BIT_OFFSET((i)))
#define GET_BIT(bm, i) (bm)[WORD_OFFSET((i))] & (((uint64_t) 1)<<BIT_OFFSET((i)))

static void
fill_bitmap_from_queue (uint64_t * restrict bm, int64_t * restrict vlist,
			int64_t out, int64_t in)
{
  for (int64_t q_index = out; q_index < in; q_index++)
    SET_BIT(bm, vlist[q_index]);
}

static void
fill_queue_from_bitmap (uint64_t * restrict bm,
			int64_t * restrict vlist,
			int64_t * restrict out, int64_t * restrict in)
{
  int64_t len = 0;

  *out = 0;

  for (int64_t k = 0; k < nv; k += 64) {
    const uint64_t m = bm[k/64];
    if (!m) continue;
    for (int kk = 0; kk < 64; ++kk)
      if (m & (((uint64_t)1) << kk))
	vlist[len++] = k + kk;
  }

  *in = len;
}

static int64_t
bfs_bottom_up_step (int64_t * restrict bfs_tree,
		    int64_t * restrict bfs_tree_depth,
		    const uint64_t * restrict past, uint64_t * restrict next,
		    int64_t level)
{
  int64_t awake_count = 0;
  const int64_t NV = nv;

  for (int64_t i = 0; i < NV; i++) {
    if (bfs_tree[i] == -1) {
      for (int64_t vo = xoff[i]; vo < xoff[1+i]; vo++) {
	const int64_t j = xadj[vo];
	if (GET_BIT(past, j)) {
	  bfs_tree[i] = j;
	  if (bfs_tree_depth) bfs_tree_depth[i] = level;
	  SET_BIT(next, i);
	  awake_count++;
	  break;
	}
      }
    }
  }

  return awake_count;
}

static void
bfs_top_down_step (int64_t * restrict bfs_tree,
		   int64_t * restrict bfs_tree_depth,
		   int64_t * restrict vlist,
		   int64_t * head_in, int64_t * tail_in,
		   int64_t level)
{
  const int64_t tail = *tail_in;
  int64_t new_tail = tail;

  for (int64_t k = *head_in; k < tail; ++k) {
    const int64_t v = vlist[k];
    const int64_t veo = xoff[1+v];
    int64_t vo;
    for (vo = xoff[v]; vo < veo; ++vo) {
      const int64_t j = xadj[vo];
      if (bfs_tree[j] == -1) {
	assert (new_tail < nv);
	vlist[new_tail++] = j;
	bfs_tree[j] = v;
	if (bfs_tree_depth) bfs_tree_depth[j] = level;
      }
    }
  }

  *head_in = tail;
  *tail_in = new_tail;

  return;
}

int
make_bfs_tree (int64_t *bfs_tree_out, int64_t * bfs_tree_depth_out, int64_t srcvtx)
{
  int64_t * restrict bfs_tree = bfs_tree_out;
  int64_t * restrict bfs_tree_depth = bfs_tree_depth_out;
  uint64_t * bm_storage = NULL;
  const int64_t NV = nv;
  const size_t bmlen = (NV + 63) & ~63;
  uint64_t * restrict past = NULL;
  uint64_t * restrict next = NULL;
  int64_t level = 0;
  int err = 0;

  int64_t * restrict vlist = NULL;
  int64_t k1, k2;

  int64_t down_cutoff = NV / BETA;
  int64_t scout_count = xoff[1+srcvtx] - xoff[srcvtx];
  int64_t awake_count = 1;
  int64_t edges_to_check = xoff[NV];

  /* Size sanity checks. */
#if 8 != CHAR_BIT
#error "Hard-coded to support eight-bit bytes."
#endif
  assert (8 == sizeof (uint64_t));

  vlist = xmalloc_large (NV * sizeof (*vlist));
  if (!vlist) return -1;

  bm_storage = xmalloc_large (2 * bmlen * sizeof (*bm_storage));
  if (!bm_storage) {
    free (vlist);
    return -2;
  }
  past = bm_storage;
  next = &bm_storage[bmlen];

  vlist[0] = srcvtx;
  k1 = 0; k2 = 1;

  for (int64_t k = 0; k < NV; ++k) {
    bfs_tree[k] = -1;
    if (bfs_tree_depth) bfs_tree_depth[k] = -1;
  }
  bfs_tree[srcvtx] = srcvtx;
  if (bfs_tree_depth) bfs_tree_depth[srcvtx] = level;

  while (awake_count != 0) {
    if (scout_count < ((edges_to_check - scout_count)/ALPHA)) {
      // Top-down
      ++level;
      bfs_top_down_step (bfs_tree, bfs_tree_depth, vlist, &k1, &k2, level);
      edges_to_check -= scout_count;
      awake_count = k2-k1;
    } else {
      // Bottom-up
      for (int64_t k = 0; k < bmlen; ++k)
	next[k] = 0;
      fill_bitmap_from_queue (next, vlist, k1, k2);
      do {
	/* Swap past & next */
	uint64_t * restrict t = past;
	past = next;
	next = t;
	for (int64_t k = 0; k < bmlen; ++k)
	  next[k] = 0;
	++level;
	awake_count = bfs_bottom_up_step (bfs_tree, bfs_tree_depth,
					  past, next, level);
      } while ((awake_count > down_cutoff));
      fill_queue_from_bitmap (next, vlist, &k1, &k2);
    }
    // Count the number of edges in the frontier
    scout_count = 0;
    for (int64_t i=k1; i<k2; i++) {
      int64_t v = vlist[i];
      scout_count += xoff[1+v] - xoff[v];
    }
  }

  xfree_large (bm_storage);
  xfree_large (vlist);

  return err;
}


struct four_heap {
  size_t n, nstore;
  int64_t keymax;
  int64_t * restrict key;
  int64_t * restrict loc;
  int16_t * restrict priority;
};

static inline size_t
parent (size_t k)
{
  return (k-1)/4;
}

static inline size_t
child (size_t k, size_t off)
{
  return 4*k+off+1;
}

#define CONCAT(a__,b__) a__ ## b__ ## __

#if defined(__GNUC__)
#define MAYBEUNUSED __attribute__ ((unused))
#else
#define MAYBEUNUSED
#endif

#define DEF_HP(hp__)							\
  const size_t N MAYBEUNUSED = (hp__).n;				\
  int64_t * restrict CONCAT(hp__, k) MAYBEUNUSED = (hp__).key;		\
  int64_t * restrict CONCAT(hp__, l) MAYBEUNUSED = (hp__).loc;		\
  int16_t * restrict CONCAT(hp__, p) MAYBEUNUSED = (hp__).priority;
#define KEY(hp__, k__) ((CONCAT(hp__, k))[(k__)])
#define LOC(hp__, k__) ((CONCAT(hp__, l))[(k__)])
#define PRIO(hp__, k__) ((CONCAT(hp__, p))[(k__)])

#define SWAPK(hp__, k1__, k2__)					\
  do {								\
  const int64_t newk2 = KEY(hp__, k1__);			\
  const int16_t tprio = PRIO(hp__, k1__);			\
  const int64_t newk1 = KEY(hp__, k1__) = KEY(hp__, k2__);	\
  PRIO(hp__, k1__) = PRIO(hp__, k2__);				\
  KEY(hp__, k2__) = newk2;					\
  PRIO(hp__, k2__) = tprio;					\
  LOC(hp__, newk1) = k1__;					\
  LOC(hp__, newk2) = k2__;					\
} while (0)

#if !defined(NDEBUG)
static inline int
check_heap (const struct four_heap hp)
{
  int out = 1;
  DEF_HP(hp);
  const int64_t KM = hp.keymax;
  for (int64_t k = 0; k < N; ++k)
    for (int ko = 0; ko < 4; ++ko) {
      int64_t kk = child (k, ko);
      if (kk < N)
	assert (PRIO(hp, k) <= PRIO(hp, kk));
    }
  for (int64_t k = 0; k < KM; ++k) {
    assert (LOC(hp, k) < N);
    assert (LOC(hp, k) < 0 || KEY(hp, LOC(hp, k)) == k);
  }
  return out;
}
#endif

static inline void
sift_up (struct four_heap hp, size_t k)
{
  DEF_HP(hp);

  while (k > 0) {
    size_t prnt = parent (k);
    if (PRIO(hp, prnt) < PRIO(hp, k)) break;
    SWAPK(hp, prnt, k);
    k = prnt;
  }
  assert (check_heap (hp));
}

static inline void
sift_down (struct four_heap hp, size_t k)
{
  DEF_HP(hp);
  int64_t childk = child (k, 0);

  while (childk < N) {
    int64_t leastk = k;
    int16_t leastprio = PRIO(hp, k);

    for (int kk = 0; kk < 4 && childk+kk < N; ++kk)
      if (PRIO(hp, childk+kk) < leastprio) {
	leastk = childk+kk;
	leastprio = PRIO(hp, leastk);
      }

    if (k == leastk) break;

    assert (leastk < N);
    SWAPK(hp, k, leastk);
    k = leastk;
    childk = child (k, 0);
  }
  assert (check_heap (hp));
}

static inline void
alloc_four_heap (struct four_heap *hp, size_t new_nstore, int64_t keymax)
{
  int64_t * new_k = xrealloc (hp->key, new_nstore * sizeof (*new_k));
  int64_t * new_loc = hp->loc;
  int16_t * new_p = xrealloc (hp->priority, new_nstore * sizeof (*new_p));
  hp->key = new_k;
  if (hp->keymax != keymax) {
    new_loc = xrealloc (hp->loc, keymax * sizeof (*new_k));
    for (int64_t k = hp->keymax; k < keymax; ++k)
      new_loc[k] = -1;
  }
  hp->loc = new_loc;
  hp->priority = new_p;
  hp->nstore = new_nstore;
  hp->keymax = keymax;
}

static inline void
init_four_heap (struct four_heap *hp, size_t nstore, int64_t keymax)
{
  memset (hp, 0, sizeof (*hp));
  alloc_four_heap (hp, nstore, keymax);
}

static inline void
fini_four_heap (struct four_heap *hp)
{
  if (!hp) return;
  if (hp->key) free (hp->key);
  if (hp->loc) free (hp->loc);
  if (hp->priority) free (hp->priority);
  memset (hp, 0, sizeof (*hp));
}

static inline void
insert (struct four_heap *hp, int64_t key, int16_t prio)
{
  if (hp->n == hp->nstore) alloc_four_heap (hp, 2*hp->nstore, hp->keymax);
  size_t where = hp->n++;
  hp->key[where] = key;
  hp->priority[where] = prio;
  hp->loc[key] = where;
  sift_up (*hp, where);
}

static inline void
find_and_decrease_prio (struct four_heap hp, int64_t key, int16_t new_prio)
{
  DEF_HP(hp);
  int64_t k = LOC(hp, key);

  if (k < 0) return;

  PRIO(hp, k) = new_prio;
  sift_up (hp, k);
}

static inline int64_t
pop_head (struct four_heap *hp)
{
  size_t n = hp->n;
  int64_t out;
  if (n == 0) return -1;
  out = hp->key[0];
  hp->loc[out] = -1;
  if (n == 1) { hp->n = 0; return out; }
#if !defined(NDEBUG)
  for (int k = 0; k < 4; ++k)
    assert (1+k >= hp->n || hp->priority[0] <= hp->priority[1+k]);
#endif
  hp->key[0] = hp->key[n-1];
  hp->priority[0] = hp->priority[n-1];
  hp->loc[hp->key[0]] = 0;
  --hp->n;
  sift_down (*hp, 0);
  return out;
}

static inline void
sssp_scan_vtx (const int64_t v, const int16_t * restrict caliber,
	       int64_t * restrict tree, int64_t * restrict dist,
	       struct four_heap * restrict queue,
	       int64_t * restrict nfixed_in, int64_t * restrict fixed,
	       char * restrict is_fixed,
	       const int64_t mu)
{
  const int64_t vend = xoff[v+1];
  int64_t nfixed = *nfixed_in;

  assert (tree[v] >= 0);
  assert (dist[v] < INT64_MAX);

  for (int64_t vo = xoff[v]; vo < vend; vo++) {
    const int64_t j = xadj[vo];
    if (is_fixed[j]) {
      assert (dist[j] <= xval[vo] + dist[v]);
      continue;
    }

    const int64_t w = xval[vo];

    if (mu + caliber[j] >= dist[j]) {
      assert (dist[j] < INT64_MAX);
      assert (dist[j] <= w + dist[v]);
      /* Then j is fixed. */
      assert (!is_fixed[j]);
      fixed[nfixed++] = j;
      is_fixed[j] = 1;
    } else {
      /* decrease key... */
      const int64_t new_d = dist[v] + w;
      assert (!is_fixed[j]);
      if (new_d < dist[j]) {
	if (dist[j] == INT64_MAX) {
	  /* enqueue */
	  insert (queue, j, new_d);
	  assert (check_heap (*queue));
	} else {
	  /* adjust */
	  find_and_decrease_prio (*queue, j, new_d);
	  assert (check_heap (*queue));
	}
	dist[j] = new_d;
	tree[j] = v;
      }
    }
  }
  *nfixed_in = nfixed;
}

int64_t max_qlen;

int
make_sssp_tree (int64_t *sssp_tree_out, int64_t * sssp_tree_dist_out,
		int64_t srcvtx)
{
  int64_t * restrict tree = sssp_tree_out;
  int64_t * restrict dist = sssp_tree_dist_out;
  const int64_t NV = nv;

  int64_t nfixed, mu;
  struct four_heap queue;
  int64_t * restrict fixed = 0;
  int16_t * restrict caliber;
  char * restrict is_fixed;

  max_qlen = 0;

  init_four_heap (&queue, NV, NV);
  fixed = xmalloc (NV * sizeof (*fixed) +
		   NV * (sizeof (*caliber) + sizeof (*is_fixed)));
  caliber = (int16_t*)&fixed[NV];
  is_fixed = (char*)&caliber[NV];

  for (int64_t k = 0; k < NV; ++k) {
    tree[k] = -1;
    dist[k] = INT64_MAX;
    is_fixed[k] = 0;
    caliber[k] = INT16_MAX;
  }

  for (int64_t k = 0; k < NV; ++k) {
    const int64_t kkend = xoff[k+1];
    int16_t c = caliber[k];
    for (int64_t kk = xoff[k]; kk < kkend; ++kk)
      if (xval[kk] < c) c = xval[kk];
    caliber[k] = c;
  }

  //fixed[0] = srcvtx;
  is_fixed[srcvtx] = 1;
  nfixed = 0;

  tree[srcvtx] = srcvtx;
  dist[srcvtx] = 0;
  mu = 0;

  sssp_scan_vtx (srcvtx, caliber, tree, dist, &queue,
		 &nfixed, fixed, is_fixed, mu);

  while (nfixed || queue.n) {
    if (queue.n > max_qlen) max_qlen = queue.n;
    if (nfixed) {
      int64_t nfixed_new, nfixed_start = 0;
      do {
	nfixed_new = nfixed;
	for (int64_t k = nfixed_start; k < nfixed; ++k) {
	  int64_t v = fixed[k];
	  sssp_scan_vtx (v, caliber, tree, dist, &queue,
			 &nfixed_new, fixed, is_fixed, mu);
	}
	nfixed_start = nfixed;
	nfixed = nfixed_new;

      } while (nfixed_new != nfixed_start); /* Exhaust the fixed vertices. */
      nfixed = 0;
    }

    assert (!nfixed);

    /* No fixed vertices, so dequeue someone. */
    if (queue.n) {
      int64_t minvtx = -1;

      do
	minvtx = pop_head (&queue);
      while (minvtx >= 0 && is_fixed[minvtx]);

      if (minvtx == -1) break;

      mu = dist[minvtx];

      sssp_scan_vtx (minvtx, caliber, tree, dist, &queue,
		     &nfixed, fixed, is_fixed, mu);
    }
  }

#if !defined(NDEBUG)
  int64_t nbad = 0;
  for (int64_t k = 0; k < NV; ++k) {
    nbad += dist[k] >= INT64_MAX;
  }
  assert (!nbad);
#endif

  free (fixed);
  fini_four_heap (&queue);

  return 0;
}

void
destroy_graph (void)
{
  free_graph ();
}
