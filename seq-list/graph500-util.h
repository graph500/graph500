/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#if !defined(GRAPH500_UTIL_HEADER_)
#define GRAPH500_UTIL_HEADER_

/** Initialize clocks. */
void init_tictoc (void);

/** Start timing. */
void tic (void);

/** Return seconds since last tic. */
double toc (void);

/** Macro to time a block. */
#define TIME(timevar, what) do { tic (); what; timevar = toc(); } while (0)

/** Compute quantiles, mean, and std. deviation. */
void statistics (double *, double*, int64_t);

/** Output timings and statistics. */
void output_results (const int64_t SCALE, int64_t nvtx_scale, int64_t edgefactor,
		     const double A, const double B, const double C, const double D,
		     const double construction_time,
		     const int NBFS, const double *bfs_time, const int64_t *bfs_nedge);

/** Allocate, aborting if the allocation fails. */
void *xmalloc (size_t);

/** Allocate a large region, aborting if the allocation fails. */
void *xmalloc_large (size_t);

/** Free memory allocated by xmalloc_large. */
void xfree_large (void *, size_t);

/** Allocate a *single* large region with backing store, aborting if the allocation fails. */
void *xmalloc_large_ext (size_t);

/** Free memory allocated by xmalloc_large_ext. */
void xfree_large_ext (void);

/** Mark external segment as unused for a bit. */
void mark_ext_unused (void);

/** Mark external segment as soon to be used. */
void mark_ext_willuse (void);

#endif /* GRAPH500_UTIL_HEADER_ */
