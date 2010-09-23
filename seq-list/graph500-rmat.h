/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#if !defined(GRAPH500_RMAT_HEADER_)
#define GRAPH500_RMAT_HEADER_

/** Fill IJ with a randomly permuted R-MAT generated edge list. */
void rmat_edgelist (int64_t *IJ, int64_t nedge, int SCALE,
		    double A, double B, double C, double D);

#endif /* GRAPH500_RMAT_HEADER_ */
