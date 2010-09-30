/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>
#include <stdio.h>

#include <mpi.h>
#include "make_graph.h"

int main(int argc, char* argv[]) {
  int log_numverts;
  int size, rank;
  unsigned long my_edges;
  unsigned long global_edges;
  double start, stop;
  size_t i;

  MPI_Init(&argc, &argv);

  log_numverts = 16; /* In base GRAPHGEN_INITIATOR_SIZE */
  if (argc >= 2) log_numverts = atoi(argv[1]);

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) fprintf(stderr, "Graph size is %" PRIu64 " vertices and %" PRIu64 " edges\n", (uint64_t)pow(2., log_numverts), (uint64_t)(pow(2., log_numverts) * 8.));

  /* Start of graph generation timing */
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  int64_t nedges;
  int64_t* result;
  double initiator[] = {.57, .19, .19, .05};
  make_graph(log_numverts, 8. * pow(2., log_numverts), 1, 2, initiator, &nedges, &result);
  MPI_Barrier(MPI_COMM_WORLD);
  stop = MPI_Wtime();
  /* End of graph generation timing */

  my_edges = 0;
  for (i = 0; i < nedges; ++i) {
    if (result[2 * i] != (uint64_t)(-1)) {
      ++my_edges;
    }
  }
  
  free(result);

  MPI_Reduce(&my_edges, &global_edges, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    fprintf(stderr, "%lu edge%s generated and permuted in %fs (%f Medges/s on %d processor%s)\n", global_edges, (global_edges == 1 ? "" : "s"), (stop - start), global_edges / (stop - start) * 1.e-6, size, (size == 1 ? "" : "s"));
  }
  MPI_Finalize();
  return 0;
}
