/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stdlib.h>

#include "xalloc.h"
#include "prng.h"
#include "generator/splittable_mrg.h"
#include "generator/permutation_gen.h"
#include "generator/graph_generator.h"

void
kronecker_edgelist (int64_t *IJ_in, int64_t nedge, int64_t SCALE,
		    double A, double B, double C)
{
  const int64_t nvtx = 1L<<SCALE;
  mrg_state state;
  int64_t * restrict IJ = IJ_in;
  int64_t * restrict vperm = NULL;
  int64_t k;
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

  rand_sort_shared ((mrg_state*)prng_state, nvtx, vperm);

  OMP("omp parallel for")
    for (k = 0; k < 2*nedge; ++k)
      if (IJ[k] >= 0)
	IJ[k] = vperm[IJ[k]];

  xfree_large (vperm);
}
