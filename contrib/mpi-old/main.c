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

#include "../graph500.h"
#include "../globals.h"
#include "../generator.h"
#include "../prng.h"
#include "../output_results.h"
#include "../xalloc.h"

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  setup_globals();
  init_prng ();

  /* Parse arguments. */
  int SCALE = 16;
  int edgefactor = 16; /* nedges / nvertices, i.e., 2*avg. degree */
  if (argc >= 2) {
    if (!strcmp (argv[1], "-v")) {
      extern char IMPLEMENTATION[];
      printf ("Graph500 (%s), code version %s, specification version %s\n",
	      IMPLEMENTATION,
	      CODE_VERSION,
	      SPEC_VERSION);
      MPI_Finalize ();
      exit (EXIT_SUCCESS);
    }
    SCALE = atoi(argv[1]);
  }
  if (argc >= 3) edgefactor = atoi(argv[2]);
  if (argc <= 1 || argc >= 4 || SCALE == 0 || edgefactor == 0) {
    if (rank == 0) {
      fprintf(stderr, "Usage: %s SCALE edgefactor\n  SCALE = log_2(# vertices) [integer, required]\n  edgefactor = (# edges) / (# vertices) = .5 * (average vertex degree) [integer, defaults to 16]\n(Random number seed and Kronecker initiator are in main.c)\n", argv[0]);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  const char* filename = getenv("TMPFILE");
  /* If filename is NULL, store data in memory */

  tuple_graph tg;
  tg.nglobaledges = (int64_t)(edgefactor) << SCALE;
  int64_t nglobalverts = (int64_t)(1) << SCALE;

  init_globals (SCALE, edgefactor, MAXWEIGHT_DEFAULT, NROOT_DEFAULT,
		A_DEFAULT, B_DEFAULT, NOISEFACT_DEFAULT);

  tg.data_in_file = (filename != NULL);

  if (tg.data_in_file) {
    MPI_File_set_errhandler(MPI_FILE_NULL, MPI_ERRORS_ARE_FATAL);
    MPI_File_open(MPI_COMM_WORLD, (char*)filename, MPI_MODE_RDWR | MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &tg.edgefile);
    MPI_File_set_size(tg.edgefile, tg.nglobaledges * sizeof(packed_edge));
    MPI_File_set_view(tg.edgefile, 0, packed_edge_mpi_type, packed_edge_mpi_type, "native", MPI_INFO_NULL);
    MPI_File_set_atomicity(tg.edgefile, 0);
  }

  /* Make the raw graph edges. */
  /* Get roots for BFS runs, plus maximum vertex with non-zero degree (used by
   * validator). */
  int64_t* bfs_roots = (int64_t*)xmalloc(NROOT * sizeof(int64_t));
  int64_t max_used_vertex = 0;

  double make_graph_start = MPI_Wtime();
  {
    /* As the graph is being generated, also keep a bitmap of vertices with
     * incident edges.  We keep a grid of processes, each row of which has a
     * separate copy of the bitmap (distributed among the processes in the
     * row), and then do an allreduce at the end.  This scheme is used to avoid
     * non-local communication and reading the file separately just to find BFS
     * roots. */
    MPI_Offset nchunks_in_file = (tg.nglobaledges + FILE_CHUNKSIZE - 1) / FILE_CHUNKSIZE;
    int64_t bitmap_size_in_bytes = int64_min(BITMAPSIZE, (nglobalverts + CHAR_BIT - 1) / CHAR_BIT);
    if (bitmap_size_in_bytes * size * CHAR_BIT < nglobalverts) {
      bitmap_size_in_bytes = (nglobalverts + size * CHAR_BIT - 1) / (size * CHAR_BIT);
    }
    int ranks_per_row = ((nglobalverts + CHAR_BIT - 1) / CHAR_BIT + bitmap_size_in_bytes - 1) / bitmap_size_in_bytes;
    int nrows = size / ranks_per_row;
    int my_row = -1, my_col = -1;
    MPI_Comm cart_comm;
    {
      int dims[2] = {size / ranks_per_row, ranks_per_row};
      int periods[2] = {0, 0};
      MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);
    }
    if (cart_comm != MPI_COMM_NULL) {
      {
        int dims[2], periods[2], coords[2];
        MPI_Cart_get(cart_comm, 2, dims, periods, coords);
        my_row = coords[0];
        my_col = coords[1];
      }
      MPI_Comm this_col;
      MPI_Comm_split(cart_comm, my_col, my_row, &this_col);
      MPI_Comm_free(&cart_comm);
      /* Every rank in a given row creates the same vertices (for updating the
       * bitmap); only one writes them to the file (or final memory buffer). */
      packed_edge* buf = (packed_edge*)xmalloc(FILE_CHUNKSIZE * sizeof(packed_edge));
      MPI_Offset block_limit = (nchunks_in_file + nrows - 1) / nrows;
      /* fprintf(stderr, "%d: nchunks_in_file = %" PRId64 ", block_limit = %" PRId64 " in grid of %d rows, %d cols\n", rank, (int64_t)nchunks_in_file, (int64_t)block_limit, nrows, ranks_per_row); */
      if (tg.data_in_file) {
        tg.edgememory_size = 0;
        tg.edgememory = NULL;
      } else {
        int my_pos = my_row + my_col * nrows;
        int last_pos = (tg.nglobaledges % ((int64_t)FILE_CHUNKSIZE * nrows * ranks_per_row) != 0) ?
                       (tg.nglobaledges / FILE_CHUNKSIZE) % (nrows * ranks_per_row) :
                       -1;
        int64_t edges_left = tg.nglobaledges % FILE_CHUNKSIZE;
        int64_t nedges = FILE_CHUNKSIZE * (tg.nglobaledges / ((int64_t)FILE_CHUNKSIZE * nrows * ranks_per_row)) +
                         FILE_CHUNKSIZE * (my_pos < (tg.nglobaledges / FILE_CHUNKSIZE) % (nrows * ranks_per_row)) +
                         (my_pos == last_pos ? edges_left : 0);
        /* fprintf(stderr, "%d: nedges = %" PRId64 " of %" PRId64 "\n", rank, (int64_t)nedges, (int64_t)tg.nglobaledges); */
        tg.edgememory_size = nedges;
        tg.edgememory = (packed_edge*)xmalloc(nedges * sizeof(packed_edge));
      }
      MPI_Offset block_idx;
      for (block_idx = 0; block_idx < block_limit; ++block_idx) {
        /* fprintf(stderr, "%d: On block %d of %d\n", rank, (int)block_idx, (int)block_limit); */
        MPI_Offset start_edge_index = int64_min(FILE_CHUNKSIZE * (block_idx * nrows + my_row), tg.nglobaledges);
        MPI_Offset edge_count = int64_min(tg.nglobaledges - start_edge_index, FILE_CHUNKSIZE);
        packed_edge* actual_buf = (!tg.data_in_file && block_idx % ranks_per_row == my_col) ?
                                  tg.edgememory + FILE_CHUNKSIZE * (block_idx / ranks_per_row) :
                                  buf;
        /* fprintf(stderr, "%d: My range is [%" PRId64 ", %" PRId64 ") %swriting into index %" PRId64 "\n", rank, (int64_t)start_edge_index, (int64_t)(start_edge_index + edge_count), (my_col == (block_idx % ranks_per_row)) ? "" : "not ", (int64_t)(FILE_CHUNKSIZE * (block_idx / ranks_per_row))); */
        if (!tg.data_in_file && block_idx % ranks_per_row == my_col) {
          assert (FILE_CHUNKSIZE * (block_idx / ranks_per_row) + edge_count <= tg.edgememory_size);
        }
	packed_edge_list (actual_buf, 0, start_edge_index, edge_count);
        if (tg.data_in_file && my_col == (block_idx % ranks_per_row)) { /* Try to spread writes among ranks */
          MPI_File_write_at(tg.edgefile, start_edge_index, actual_buf, edge_count, packed_edge_mpi_type, MPI_STATUS_IGNORE);
        }
      }
      free(buf);
      MPI_Comm_free(&this_col);
    } else {
      tg.edgememory = NULL;
      tg.edgememory_size = 0;
    }
    MPI_Allreduce(&tg.edgememory_size, &tg.max_edgememory_size, 1, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);
    /* Find roots and max used vertex */
    sample_roots (bfs_roots, NROOT, NE);
    if (tg.data_in_file) {
      MPI_File_sync(tg.edgefile);
    }
  }
  max_used_vertex = NV-1;
  double make_graph_stop = MPI_Wtime();
  double make_graph_time = make_graph_stop - make_graph_start;
  if (rank == 0) { /* Not an official part of the results */
    if (getenv("VERBOSE"))
      fprintf(stderr, "graph_generation:               %f s\n", make_graph_time);
  }

  /* Make user's graph data structure. */
  double data_struct_start = MPI_Wtime();
  make_graph_data_structure(&tg);
  double data_struct_stop = MPI_Wtime();
  double data_struct_time = data_struct_stop - data_struct_start;
  if (rank == 0) { /* Not an official part of the results */
    if (getenv("VERBOSE"))
      fprintf(stderr, "construction_time:              %f s\n", data_struct_time);
  }

  /* Number of edges visited in each BFS; a double so get_statistics can be
   * used directly. */
  double* edge_counts = (double*)xmalloc(NROOT * sizeof(double));

  /* Run BFS. */
  int validation_passed = 1;
  double* bfs_times = (double*)xmalloc(NROOT * sizeof(double));
  int64_t* bfs_depth = (int64_t*)xmalloc(NROOT * sizeof(*bfs_depth));
  double* validate_times = (double*)xmalloc(NROOT * sizeof(double));
  uint64_t nlocalverts = get_nlocalverts_for_pred();
  int64_t* pred = (int64_t*)xMPI_Alloc_mem(nlocalverts * sizeof(int64_t));

  int bfs_root_idx;
  for (bfs_root_idx = 0; bfs_root_idx < NROOT; ++bfs_root_idx) {
    int64_t root = bfs_roots[bfs_root_idx];

    if (rank == 0 && getenv ("VERBOSE"))
      fprintf(stderr, "Running BFS %d\n", bfs_root_idx);

    /* Clear the pred array. */
    memset(pred, 0, nlocalverts * sizeof(int64_t));

    /* Do the actual BFS. */
    double bfs_start = MPI_Wtime();
    run_bfs(root, &pred[0]);
    double bfs_stop = MPI_Wtime();
    bfs_times[bfs_root_idx] = bfs_stop - bfs_start;
    if (rank == 0 && getenv ("VERBOSE"))
      fprintf(stderr, "Time for BFS %d is %f\n", bfs_root_idx, bfs_times[bfs_root_idx]);

    /* Validate result. */
    if (rank == 0 && getenv ("VERBOSE"))
      fprintf(stderr, "Validating BFS %d\n", bfs_root_idx);

    double validate_start = MPI_Wtime();
    int64_t edge_visit_count;
    int validation_passed_one = validate_bfs_result(&tg, max_used_vertex + 1, nlocalverts, root, pred, &edge_visit_count, &bfs_depth[bfs_root_idx]);
    double validate_stop = MPI_Wtime();
    validate_times[bfs_root_idx] = validate_stop - validate_start;
    if (rank == 0 && getenv ("VERBOSE"))
      fprintf(stderr, "Validate time for BFS %d is %f\n", bfs_root_idx, validate_times[bfs_root_idx]);
    edge_counts[bfs_root_idx] = (double)edge_visit_count;
    if (rank == 0 && getenv ("VERBOSE"))
      fprintf(stderr, "TEPS for BFS %d is %g\n", bfs_root_idx, edge_visit_count / bfs_times[bfs_root_idx]);

    if (!validation_passed_one) {
      validation_passed = 0;
      if (rank == 0 && getenv ("VERBOSE"))
	fprintf(stderr, "Validation failed for this BFS root; skipping rest.\n");
      break;
    }
  }

  MPI_Free_mem(pred);
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
      extern char IMPLEMENTATION[];

      output_results (IMPLEMENTATION,
		      make_graph_time, data_struct_time,
		      bfs_roots,
		      bfs_times, bfs_depth, validate_times);
    }
  }
  free(bfs_roots);
  free(bfs_times);
  free(validate_times);

  cleanup_globals();
  MPI_Finalize();
  return 0;
}
