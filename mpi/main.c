/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

/* These need to be before any possible inclusions of stdint.h or inttypes.h.
 * */
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#include "../generator/make_graph.h"
#include "../generator/utils.h"
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
  int edgefactor = 16; /* nedges / nvertices, i.e., 2*avg. degree */
  if (argc >= 2) SCALE = atoi(argv[1]);
  if (argc >= 3) edgefactor = atoi(argv[2]);
  if (argc <= 1 || argc >= 4 || SCALE == 0 || edgefactor == 0) {
    if (rank == 0) {
      fprintf(stderr, "Usage: %s SCALE edgefactor\n  SCALE = log_2(# vertices) [integer, required]\n  edgefactor = (# edges) / (# vertices) = .5 * (average vertex degree) [integer, defaults to 16]\n(Random number seed and Kronecker initiator are in main.c)\n", argv[0]);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  uint64_t seed1 = 2, seed2 = 3;

  const char* filename = getenv("TMPFILE");
  /* If filename is NULL, store data in memory */

  tuple_graph tg;
  tg.nglobaledges = (int64_t)(edgefactor) << SCALE;

  if (filename != NULL) {
    MPI_File_set_errhandler(MPI_FILE_NULL, MPI_ERRORS_ARE_FATAL);
    MPI_File_open(MPI_COMM_WORLD, (char*)filename, MPI_MODE_RDWR | MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &tg.edgefile);
    MPI_File_set_size(tg.edgefile, tg.nglobaledges * sizeof(packed_edge));
    MPI_File_set_view(tg.edgefile, 0, packed_edge_mpi_type, packed_edge_mpi_type, "native", MPI_INFO_NULL);
    MPI_File_set_atomicity(tg.edgefile, 0);
  }

  /* Make the raw graph edges. */
  double make_graph_start = MPI_Wtime();
  {
    /* Spread the two 64-bit numbers into five nonzero values in the correct
     * range. */
    uint_fast32_t seed[5];
    make_mrg_seed(seed1, seed2, seed);

    if (filename != NULL) {
      tg.data_in_file = 1;
      MPI_Offset nchunks_in_file = (tg.nglobaledges + FILE_CHUNKSIZE - 1) / FILE_CHUNKSIZE;
      MPI_Offset block_limit = DIV_SIZE(nchunks_in_file + size - 1);
      MPI_Offset block_idx;
      packed_edge* buf = xmalloc(FILE_CHUNKSIZE * sizeof(packed_edge));
      for (block_idx = 0; block_idx < block_limit; ++block_idx) {
        MPI_Offset start_edge_index = FILE_CHUNKSIZE * (MUL_SIZE(block_idx) + rank);
        if (start_edge_index > tg.nglobaledges) start_edge_index = tg.nglobaledges;
        MPI_Offset edge_count = tg.nglobaledges - start_edge_index;
        if (edge_count > FILE_CHUNKSIZE) edge_count = FILE_CHUNKSIZE;
        generate_kronecker_range(seed, SCALE, start_edge_index, start_edge_index + edge_count, buf);
        MPI_File_write_at_all(tg.edgefile, start_edge_index, buf, edge_count, packed_edge_mpi_type, MPI_STATUS_IGNORE);
      }
      free(buf);
      MPI_File_sync(tg.edgefile);
      tg.edgememory_size = 0;
      tg.edgememory = NULL;
    } else {
      tg.data_in_file = 0;
      int64_t start_idx = rank * DIV_SIZE(tg.nglobaledges) + (rank < MOD_SIZE(tg.nglobaledges) ? rank : MOD_SIZE(tg.nglobaledges));
      int64_t end_idx = (rank + 1) * DIV_SIZE(tg.nglobaledges) + (rank + 1 < MOD_SIZE(tg.nglobaledges) ? rank + 1 : MOD_SIZE(tg.nglobaledges));
      int64_t nedges = end_idx - start_idx;
      tg.edgememory_size = nedges;
      tg.edgememory = (packed_edge*)xmalloc(nedges * sizeof(packed_edge));
      generate_kronecker_range(seed, SCALE, start_idx, end_idx, tg.edgememory);
    }
  }
  double make_graph_stop = MPI_Wtime();
  double make_graph_time = make_graph_stop - make_graph_start;
  if (rank == 0) { /* Not an official part of the results */
    fprintf(stderr, "graph_generation:               %f s\n", make_graph_time);
  }

  /* Get roots for BFS runs, plus maximum vertex with non-zero degree (used by
   * validator). */
  // int num_bfs_roots = 64;
  int num_bfs_roots = 16;
  int64_t* bfs_roots = (int64_t*)xmalloc(num_bfs_roots * sizeof(int64_t));
  int64_t max_used_vertex = 0;
  double find_roots_start = MPI_Wtime();
  find_bfs_roots(&num_bfs_roots, (INT64_C(1) << SCALE), &tg, seed1, seed2, bfs_roots, &max_used_vertex);
  double find_roots_stop = MPI_Wtime();
  double find_roots_time = find_roots_stop - find_roots_start;
  if (rank == 0) { /* Not an official part of the results */
    fprintf(stderr, "find_roots:                     %f s\n", find_roots_time);
  }

  /* Make user's graph data structure. */
  double data_struct_start = MPI_Wtime();
  make_graph_data_structure(&tg);
  double data_struct_stop = MPI_Wtime();
  double data_struct_time = data_struct_stop - data_struct_start;
  if (rank == 0) { /* Not an official part of the results */
    fprintf(stderr, "construction_time:              %f s\n", data_struct_time);
  }

  /* Number of edges visited in each BFS; a double so get_statistics can be
   * used directly. */
  double* edge_counts = (double*)xmalloc(num_bfs_roots * sizeof(double));

  /* Run BFS. */
  int validation_passed = 1;
  double* bfs_times = (double*)xmalloc(num_bfs_roots * sizeof(double));
  double* validate_times = (double*)xmalloc(num_bfs_roots * sizeof(double));
  uint64_t nlocalverts = get_nlocalverts_for_pred();
  int64_t* pred = (int64_t*)xMPI_Alloc_mem(nlocalverts * sizeof(int64_t));

  int bfs_root_idx;
  for (bfs_root_idx = 0; bfs_root_idx < num_bfs_roots; ++bfs_root_idx) {
    int64_t root = bfs_roots[bfs_root_idx];

    if (rank == 0) fprintf(stderr, "Running BFS %d\n", bfs_root_idx);

    /* Clear the pred array. */
    memset(pred, 0, nlocalverts * sizeof(int64_t));

    /* Do the actual BFS. */
    double bfs_start = MPI_Wtime();
    run_bfs(root, &pred[0]);
    double bfs_stop = MPI_Wtime();
    bfs_times[bfs_root_idx] = bfs_stop - bfs_start;
    if (rank == 0) fprintf(stderr, "Time for BFS %d is %f\n", bfs_root_idx, bfs_times[bfs_root_idx]);

    /* Validate result. */
    if (rank == 0) fprintf(stderr, "Validating BFS %d\n", bfs_root_idx);

    double validate_start = MPI_Wtime();
    int64_t edge_visit_count;
    int validation_passed_one = validate_bfs_result(&tg, max_used_vertex + 1, nlocalverts, root, pred, &edge_visit_count);
    double validate_stop = MPI_Wtime();
    validate_times[bfs_root_idx] = validate_stop - validate_start;
    if (rank == 0) fprintf(stderr, "Validate time for BFS %d is %f\n", bfs_root_idx, validate_times[bfs_root_idx]);
    edge_counts[bfs_root_idx] = (double)edge_visit_count;
    if (rank == 0) fprintf(stderr, "TEPS for BFS %d is %g\n", bfs_root_idx, edge_visit_count / bfs_times[bfs_root_idx]);

    if (!validation_passed_one) {
      validation_passed = 0;
      if (rank == 0) fprintf(stderr, "Validation failed for this BFS root; skipping rest.\n");
      break;
    }
  }

  MPI_Free_mem(pred);
  free(bfs_roots);
  free_graph_data_structure();

  if (tg.data_in_file) {
    MPI_File_close(&tg.edgefile);
  } else {
    free(tg.edgememory); tg.edgememory = NULL;
  }

  /* Print results. */
  if (rank == 0) {
    if (!validation_passed) {
      fprintf(stdout, "No results printed for invalid run.\n");
    } else {
      int i;
      fprintf(stdout, "SCALE:                          %d\n", SCALE);
      fprintf(stdout, "edgefactor:                     %d\n", edgefactor);
      fprintf(stdout, "NBFS:                           %d\n", num_bfs_roots);
      fprintf(stdout, "graph_generation:               %g\n", make_graph_time);
      fprintf(stdout, "num_mpi_processes:              %d\n", size);
      fprintf(stdout, "construction_time:              %g\n", data_struct_time);
      double stats[s_LAST];
      get_statistics(bfs_times, num_bfs_roots, stats);
      fprintf(stdout, "min_time:                       %g\n", stats[s_minimum]);
      fprintf(stdout, "firstquartile_time:             %g\n", stats[s_firstquartile]);
      fprintf(stdout, "median_time:                    %g\n", stats[s_median]);
      fprintf(stdout, "thirdquartile_time:             %g\n", stats[s_thirdquartile]);
      fprintf(stdout, "max_time:                       %g\n", stats[s_maximum]);
      fprintf(stdout, "mean_time:                      %g\n", stats[s_mean]);
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
      fprintf(stdout, "min_TEPS:                       %g\n", 1. / stats[s_maximum]);
      fprintf(stdout, "firstquartile_TEPS:             %g\n", 1. / stats[s_thirdquartile]);
      fprintf(stdout, "median_TEPS:                    %g\n", 1. / stats[s_median]);
      fprintf(stdout, "thirdquartile_TEPS:             %g\n", 1. / stats[s_firstquartile]);
      fprintf(stdout, "max_TEPS:                       %g\n", 1. / stats[s_minimum]);
      fprintf(stdout, "harmonic_mean_TEPS:             %g\n", 1. / stats[s_mean]);
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
      fprintf(stdout, "min_validate:                   %g\n", stats[s_minimum]);
      fprintf(stdout, "firstquartile_validate:         %g\n", stats[s_firstquartile]);
      fprintf(stdout, "median_validate:                %g\n", stats[s_median]);
      fprintf(stdout, "thirdquartile_validate:         %g\n", stats[s_thirdquartile]);
      fprintf(stdout, "max_validate:                   %g\n", stats[s_maximum]);
      fprintf(stdout, "mean_validate:                  %g\n", stats[s_mean]);
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

  cleanup_globals();
  MPI_Finalize();
  return 0;
}
