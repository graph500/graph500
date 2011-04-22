/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stdlib.h>

#include "xalloc.h"
#include "prng.h"
#include "rmat.h"
#include "generator/splittable_mrg.h"
#include "generator/graph_generator.h"

#if 0
void
kronecker_edgelist (struct packed_edge *IJ_in, int64_t nedge, int64_t SCALE,
		    double A, double B, double C)
{
  const int64_t nvtx = 1L<<SCALE;
  struct packed_edge * restrict IJ = IJ_in;
  int64_t * restrict vperm = NULL;
  double D;

  D = 1.0 - (A+B+C);

  vperm = xmalloc_large_ext (nvtx * sizeof (*vperm));

#if !defined(__MTA__)
  double initiator[] = {A, B, C, D};

  OMP("omp parallel") {
    generate_kronecker (omp_get_thread_num (), omp_get_num_threads (),
			prng_seed, SCALE, nedge, initiator, IJ);
  }
#else /* __MTA__ */
  int64_t rank, size;
  MTA("mta use 100 streams") MTA("mta for all streams rank of size") {
    double tinitiator[] = {A, B, C, D};
    generate_kronecker (rank, size, prng_seed, SCALE, nedge, tinitiator, IJ);
  }
#endif

  OMP("omp parallel") {
    permute_vertex_labels (IJ, nedge, nvtx, prng_state, vperm);

    OMP("omp barrier");
    OMP("omp master") xfree_large (vperm);

    permute_edgelist (IJ, nedge, prng_state);
  }
}
#endif
