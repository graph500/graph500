/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#define GRAPHGEN_DISTRIBUTED_MEMORY
#define GRAPH_GENERATOR_MPI

/* These need to be before any possible inclusions of stdint.h or inttypes.h.
 * */
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#include "../generator/make_graph.h"
#include "common.h"
#include <math.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>

static int compare_doubles(const void* a, const void* b) {
  double aa = *(const double*)a;
  double bb = *(const double*)b;
  return (aa < bb) ? -1 : (aa == bb) ? 0 : 1;
}

enum {s_minimum, s_firstquartile, s_median, s_thirdquartile, s_maximum, s_mean, s_std, s_LAST};
static void get_statistics(const double x[], int n, double r[s_LAST]) {
  double temp;
  int i;
  /* Compute mean. */
  temp = 0;
  for (i = 0; i < n; ++i) temp += x[i];
  temp /= n;
  r[s_mean] = temp;
  /* Compute std. dev. */
  temp = 0;
  for (i = 0; i < n; ++i) temp += (x[i] - r[s_mean]) * (x[i] - r[s_mean]);
  temp /= n - 1;
  r[s_std] = sqrt(temp);
  /* Sort x. */
  double* xx = (double*)xmalloc(n * sizeof(double));
  memcpy(xx, x, n * sizeof(double));
  qsort(xx, n, sizeof(double), compare_doubles);
  /* Get order statistics. */
  r[s_minimum] = xx[0];
  r[s_firstquartile] = (xx[(n - 1) / 4] + xx[n / 4]) * .5;
  r[s_median] = (xx[(n - 1) / 2] + xx[n / 2]) * .5;
  r[s_thirdquartile] = (xx[n - 1 - (n - 1) / 4] + xx[n - 1 - n / 4]) * .5;
  r[s_maximum] = xx[n - 1];
  /* Clean up. */
  free(xx);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  setup_globals();

  /* Parse arguments. */
  int SCALE = 16;
  double edgefactor = 16.; /* nedges / nvertices, i.e., 2*avg. degree */
  if (argc >= 2) SCALE = atoi(argv[1]);
  if (argc >= 3) edgefactor = atof(argv[2]);
  if (argc <= 1 || argc >= 4 || SCALE == 0 || edgefactor == 0) {
    if (rank == 0) {
      fprintf(stderr, "Usage: %s SCALE edgefactor\n  SCALE = log_2(# vertices) [integer, required]\n  edgefactor = (# edges) / (# vertices) = .5 * (average vertex degree) [float, defaults to 16]\n(Random number seed and Kronecker initiator are in main.c)\n", argv[0]);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  uint64_t seed1 = 2, seed2 = 3;
  const double initiator[4] = {.57, .19, .19, .05};

  int64_t nedges;
  int64_t* edges;

  /* Make the raw graph edges. */
  double make_graph_start = MPI_Wtime();
  make_graph(SCALE, (int64_t)(edgefactor * pow(2., SCALE)), seed1, seed2, initiator,
             &nedges, &edges);
  double make_graph_stop = MPI_Wtime();
  double make_graph_time = make_graph_stop - make_graph_start;

  /* CSR graph data structure. */
  csr_graph g;

  /* Make CSR data structure, redistributing data using VERTEX_OWNER
   * distribution. */
  double data_struct_start = MPI_Wtime();
  convert_graph_to_csr(nedges, edges, &g);
  double data_struct_stop = MPI_Wtime();
  double data_struct_time = data_struct_stop - data_struct_start;
  free(edges); edges = NULL;

  /* Get roots for BFS runs. */
  int num_bfs_roots = 64;
  int64_t* bfs_roots = (int64_t*)xmalloc(num_bfs_roots * sizeof(int64_t));
  find_bfs_roots(&num_bfs_roots, &g, seed1, seed2, bfs_roots);

  /* Number of edges visited in each BFS; a double so get_statistics can be
   * used directly. */
  double* edge_counts = (double*)xmalloc(num_bfs_roots * sizeof(double));

  /* Run BFS. */
  int validation_passed = 1;
  double* bfs_times = (double*)xmalloc(num_bfs_roots * sizeof(double));
  double* validate_times = (double*)xmalloc(num_bfs_roots * sizeof(double));
  int64_t* pred = (int64_t*)xMPI_Alloc_mem(g.nlocalverts * sizeof(int64_t));

  int bfs_root_idx;
  for (bfs_root_idx = 0; bfs_root_idx < num_bfs_roots; ++bfs_root_idx) {
    int64_t root = bfs_roots[bfs_root_idx];

    if (rank == 0) fprintf(stderr, "Running BFS %d\n", bfs_root_idx);

    /* Clear the pred array. */
    memset(pred, 0, g.nlocalverts * sizeof(int64_t));

    /* Do the actual BFS. */
    double bfs_start = MPI_Wtime();
    int64_t nvisited = 0;
    run_mpi_bfs(&g, root, &pred[0], &nvisited);
    double bfs_stop = MPI_Wtime();
    bfs_times[bfs_root_idx] = bfs_stop - bfs_start;

    /* Validate result. */
    if (rank == 0) fprintf(stderr, "Validating BFS %d\n", bfs_root_idx);

    double validate_start = MPI_Wtime();
    int validation_passed_one = validate_bfs_result(&g, root, pred, nvisited);
    double validate_stop = MPI_Wtime();
    validate_times[bfs_root_idx] = validate_stop - validate_start;

    if (!validation_passed_one) {
      validation_passed = 0;
      if (rank == 0) fprintf(stderr, "Validation failed for this BFS root; skipping rest.\n");
      break;
    }

    /* Calculate number of input edges visited. */
    {
      int64_t edge_visit_count = 0;
      size_t v_local;
      for (v_local = 0; v_local < g.nlocalverts; ++v_local) {
        if (pred[v_local] != -1) {
          size_t ei, ei_end = g.rowstarts[v_local + 1];
          for (ei = g.rowstarts[v_local]; ei < ei_end; ++ei) {
            /* Test here is so that each input edge is counted exactly once, even
             * though non-self-loops are duplicated in the CSR data structure. */
            if (g.column[ei] <= VERTEX_TO_GLOBAL(v_local)) {
              ++edge_visit_count;
            }
          }
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &edge_visit_count, 1, INT64_T_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
      edge_counts[bfs_root_idx] = (double)edge_visit_count;
    }
  }

  MPI_Free_mem(pred);
  free(bfs_roots);
  free_csr_graph(&g);

  /* Print results. */
  if (rank == 0) {
    if (!validation_passed) {
      fprintf(stdout, "No results printed for invalid run.\n");
    } else {
      int i;
      fprintf(stdout, "SCALE:                          %d\n", SCALE);
      fprintf(stdout, "edgefactor:                     %.2g\n", edgefactor);
      fprintf(stdout, "NBFS:                           %d\n", num_bfs_roots);
      fprintf(stdout, "graph_generation:               %g s\n", make_graph_time);
      fprintf(stdout, "num_mpi_processes:              %d\n", size);
      fprintf(stdout, "construction_time:              %g s\n", data_struct_time);
      double stats[s_LAST];
      get_statistics(bfs_times, num_bfs_roots, stats);
      fprintf(stdout, "min_time:                       %g s\n", stats[s_minimum]);
      fprintf(stdout, "firstquartile_time:             %g s\n", stats[s_firstquartile]);
      fprintf(stdout, "median_time:                    %g s\n", stats[s_median]);
      fprintf(stdout, "thirdquartile_time:             %g s\n", stats[s_thirdquartile]);
      fprintf(stdout, "max_time:                       %g s\n", stats[s_maximum]);
      fprintf(stdout, "mean_time:                      %g s\n", stats[s_mean]);
      fprintf(stdout, "stddev_time:                    %g\n", stats[s_std]);
      get_statistics(edge_counts, num_bfs_roots, stats);
      fprintf(stdout, "min_nedge:                      %.11g\n", stats[s_minimum]);
      fprintf(stdout, "firstquartile_nedge:            %.11g\n", stats[s_firstquartile]);
      fprintf(stdout, "median_nedge:                   %.11g\n", stats[s_median]);
      fprintf(stdout, "thirdquartile_nedge:            %.11g\n", stats[s_thirdquartile]);
      fprintf(stdout, "max_nedge:                      %.11g\n", stats[s_maximum]);
      fprintf(stdout, "mean_nedge:                     %.11g\n", stats[s_mean]);
      fprintf(stdout, "stddev_nedge:                   %.11g\n", stats[s_std]);
      double* secs_per_edge = (double*)xmalloc(num_bfs_roots * sizeof(double));
      for (i = 0; i < num_bfs_roots; ++i) secs_per_edge[i] = bfs_times[i] / edge_counts[i];
      get_statistics(secs_per_edge, num_bfs_roots, stats);
      fprintf(stdout, "min_TEPS:                       %g TEPS\n", 1. / stats[s_maximum]);
      fprintf(stdout, "firstquartile_TEPS:             %g TEPS\n", 1. / stats[s_thirdquartile]);
      fprintf(stdout, "median_TEPS:                    %g TEPS\n", 1. / stats[s_median]);
      fprintf(stdout, "thirdquartile_TEPS:             %g TEPS\n", 1. / stats[s_firstquartile]);
      fprintf(stdout, "max_TEPS:                       %g TEPS\n", 1. / stats[s_minimum]);
      fprintf(stdout, "harmonic_mean_TEPS:             %g TEPS\n", 1. / stats[s_mean]);
      /* Formula from:
       * Title: The Standard Errors of the Geometric and Harmonic Means and
       *        Their Application to Index Numbers
       * Author(s): Nilan Norris
       * Source: The Annals of Mathematical Statistics, Vol. 11, No. 4 (Dec., 1940), pp. 445-448
       * Publisher(s): Institute of Mathematical Statistics
       * Stable URL: http://www.jstor.org/stable/2235723
       * (same source as in specification). */
      fprintf(stdout, "harmonic_stddev_TEPS:           %g\n", stats[s_std] / (stats[s_mean] * stats[s_mean] * sqrt(num_bfs_roots - 1)));
      free(secs_per_edge); secs_per_edge = NULL;
      free(edge_counts); edge_counts = NULL;
      get_statistics(validate_times, num_bfs_roots, stats);
      fprintf(stdout, "min_validate:                   %g s\n", stats[s_minimum]);
      fprintf(stdout, "firstquartile_validate:         %g s\n", stats[s_firstquartile]);
      fprintf(stdout, "median_validate:                %g s\n", stats[s_median]);
      fprintf(stdout, "thirdquartile_validate:         %g s\n", stats[s_thirdquartile]);
      fprintf(stdout, "max_validate:                   %g s\n", stats[s_maximum]);
      fprintf(stdout, "mean_validate:                  %g s\n", stats[s_mean]);
      fprintf(stdout, "stddev_validate:                %g\n", stats[s_std]);
#if 0
      for (i = 0; i < num_bfs_roots; ++i) {
        fprintf(stdout, "Run %3d:                        %g s, validation %g s\n", i + 1, bfs_times[i], validate_times[i]);
      }
#endif
    }
  }
  free(bfs_times);
  free(validate_times);

  MPI_Finalize();
  return 0;
}
