/* -*- mode: C; fill-column: 70; -*- */
/* Copyright 2010-2011,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#if !defined(OUTPUT_RESULTS_HEADER_)
#define OUTPUT_RESULTS_HEADER_

void
output_results (const char * implementation,
		const double generation_time,
		const double construction_time,
		const int64_t *root,
		const double *bfs_time,
		const int64_t *bfs_depth,
		const double *bfs_verify_time,
		const double *sssp_time,
		const int64_t *sssp_depth,
		const double *sssp_verify_time);

#endif /* OUTPUT_RESULTS_HEADER_ */
