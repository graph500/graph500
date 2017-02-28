/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */
/*           Anton Korzh                                                   */

#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <stddef.h>
#include <limits.h>
#include <mpi.h>
#include "../generator/graph_generator.h"
#ifndef PROCS_PER_NODE_NOT_POWER_OF_TWO
#define SIZE_MUST_BE_A_POWER_OF_TWO
#endif
extern int rank, size;
#ifdef SIZE_MUST_BE_A_POWER_OF_TWO
extern int lgsize;
#endif
extern MPI_Datatype packed_edge_mpi_type; /* MPI datatype for packed_edge struct */

static const int ulong_bits = sizeof(unsigned long) * CHAR_BIT;

/* Distribute edges by their endpoints (make two directed copies of each input
 * undirected edge); distribution is 1-d and cyclic. */
#ifdef SIZE_MUST_BE_A_POWER_OF_TWO
#define MOD_SIZE(v) ((v) & ((1 << lgsize) - 1))
#define DIV_SIZE(v) ((v) >> lgsize)
#define MUL_SIZE(x) ((x) << lgsize)
#else
#define MOD_SIZE(v) ((v) % size)
#define DIV_SIZE(v) ((v) / size)
#define MUL_SIZE(x) ((x) * size)
#endif

#define VERTEX_OWNER(v) ((int)(MOD_SIZE(v)))
#define VERTEX_LOCAL(v) ((size_t)(DIV_SIZE(v)))
#define VERTEX_TO_GLOBAL(r, i) ((int64_t)(MUL_SIZE((uint64_t)((i))) + (int)((r))))

typedef struct tuple_graph {
	int data_in_file; /* 1 for file, 0 for memory */
	int write_file; /* 1 if the file needs written, 0 if re-used and read */
	packed_edge* restrict edgememory; /* NULL if edges are in file */
	int64_t edgememory_size;
	int64_t max_edgememory_size;
	MPI_File edgefile; /* Or MPI_FILE_NULL if edges are in memory */
	int64_t nglobaledges; /* Number of edges in graph, in both cases */
#ifdef SSSP
	float* restrict weightmemory;
	MPI_File weightfile;
#endif
} tuple_graph;

#define FILE_CHUNKSIZE ((MPI_Offset)(1) << 23) /* Size of one file I/O block or memory block to be processed in one step, in edges */

/* Simple iteration of edge data or file; cannot be nested. */
#define ITERATE_TUPLE_GRAPH_BLOCK_COUNT(tg) \
	((tg)->data_in_file ? \
	 (DIV_SIZE((MPI_Offset)(((tg)->nglobaledges + FILE_CHUNKSIZE - 1) / FILE_CHUNKSIZE) \
		   + size - 1)) : \
	 (((tg)->max_edgememory_size + FILE_CHUNKSIZE - 1) / FILE_CHUNKSIZE))
#ifdef SSSP
#define ITERATE_TUPLE_GRAPH_BEGIN(tg, user_buf, user_buf_count,wbuf) \
	do { \
		MPI_Offset block_limit = ITERATE_TUPLE_GRAPH_BLOCK_COUNT(tg); \
		if ((tg)->data_in_file) { \
			assert ((tg)->edgefile != MPI_FILE_NULL); \
			assert ((tg)->weightfile != MPI_FILE_NULL); \
		} \
		/* fprintf(stderr, "%d has block_limit = %td\n", rank, (ptrdiff_t)block_limit); */ \
		MPI_Offset block_idx; \
		packed_edge* edge_data_from_file = (packed_edge*)((tg)->data_in_file ? xmalloc(FILE_CHUNKSIZE * sizeof(packed_edge)) : NULL); \
		float* weight_data_from_file = (float*)((tg)->data_in_file ? xmalloc(FILE_CHUNKSIZE * sizeof(float)) : NULL); \
		int64_t edge_count_i = (int64_t)(-1); \
		if ((tg)->data_in_file && block_limit > 0) { \
			MPI_Offset start_edge_index = FILE_CHUNKSIZE * rank; \
			if (start_edge_index > (tg)->nglobaledges) start_edge_index = (tg)->nglobaledges; \
			edge_count_i = (tg)->nglobaledges - start_edge_index; \
			if (edge_count_i > FILE_CHUNKSIZE) edge_count_i = FILE_CHUNKSIZE; \
			MPI_File_read_at_all_begin((tg)->edgefile, start_edge_index, edge_data_from_file, edge_count_i, packed_edge_mpi_type); \
			MPI_File_read_at_all_begin((tg)->weightfile, start_edge_index, weight_data_from_file, edge_count_i, MPI_FLOAT); \
		} \
		int break_from_block_loop = 0; \
		for (block_idx = 0; block_idx < block_limit; ++block_idx) { \
			MPI_Offset start_edge_index, end_edge_index; \
			if ((tg)->data_in_file) { \
				start_edge_index = FILE_CHUNKSIZE * (MUL_SIZE(block_idx) + rank); \
				if (start_edge_index > (tg)->nglobaledges) start_edge_index = (tg)->nglobaledges; \
				end_edge_index = start_edge_index + FILE_CHUNKSIZE; \
				if (end_edge_index > (tg)->nglobaledges) end_edge_index = (tg)->nglobaledges; \
				/* fprintf(stderr, "%d trying to read offset = %" PRId64 ", count = %" PRId64 "\n", rank, start_edge_index, edge_count_i); */ \
				MPI_File_read_at_all_end((tg)->edgefile, edge_data_from_file, MPI_STATUS_IGNORE); \
				MPI_File_read_at_all_end((tg)->weightfile, weight_data_from_file, MPI_STATUS_IGNORE); \
			} else { \
				start_edge_index = int64_min(FILE_CHUNKSIZE * block_idx, (tg)->edgememory_size); \
				end_edge_index = int64_min(start_edge_index + FILE_CHUNKSIZE, (tg)->edgememory_size); \
			} \
			edge_count_i = end_edge_index - start_edge_index; \
			const packed_edge* restrict const user_buf = ((tg)->data_in_file ? edge_data_from_file : (tg)->edgememory + start_edge_index); \
			const float* restrict const wbuf = ((tg)->data_in_file ? weight_data_from_file : (tg)->weightmemory + start_edge_index); \
			ptrdiff_t const user_buf_count = edge_count_i; \
			assert (user_buf != NULL); \
			assert (wbuf != NULL); \
			assert (user_buf_count >= 0); \
			assert (tuple_graph_max_bufsize((tg)) >= user_buf_count); \
			int iteration_count = 0; (void)iteration_count; \
			int buffer_released_this_iter = 0; /* To allow explicit buffer release to be optional */ \
			while (1) { \
				/* Prevent continue */ assert (iteration_count == 0); \
				buffer_released_this_iter = 0; \
				{
#else
#define ITERATE_TUPLE_GRAPH_BEGIN(tg, user_buf, user_buf_count,unused) \
					do { \
						MPI_Offset block_limit = ITERATE_TUPLE_GRAPH_BLOCK_COUNT(tg); \
						if ((tg)->data_in_file) { \
							assert ((tg)->edgefile != MPI_FILE_NULL); \
						} \
						/* fprintf(stderr, "%d has block_limit = %td\n", rank, (ptrdiff_t)block_limit); */ \
						MPI_Offset block_idx; \
						packed_edge* edge_data_from_file = (packed_edge*)((tg)->data_in_file ? xmalloc(FILE_CHUNKSIZE * sizeof(packed_edge)) : NULL); \
						int64_t edge_count_i = (int64_t)(-1); \
						if ((tg)->data_in_file && block_limit > 0) { \
							MPI_Offset start_edge_index = FILE_CHUNKSIZE * rank; \
							if (start_edge_index > (tg)->nglobaledges) start_edge_index = (tg)->nglobaledges; \
							edge_count_i = (tg)->nglobaledges - start_edge_index; \
							if (edge_count_i > FILE_CHUNKSIZE) edge_count_i = FILE_CHUNKSIZE; \
							MPI_File_read_at_all_begin((tg)->edgefile, start_edge_index, edge_data_from_file, edge_count_i, packed_edge_mpi_type); \
						} \
						int break_from_block_loop = 0; \
						for (block_idx = 0; block_idx < block_limit; ++block_idx) { \
							MPI_Offset start_edge_index, end_edge_index; \
							if ((tg)->data_in_file) { \
								start_edge_index = FILE_CHUNKSIZE * (MUL_SIZE(block_idx) + rank); \
								if (start_edge_index > (tg)->nglobaledges) start_edge_index = (tg)->nglobaledges; \
								end_edge_index = start_edge_index + FILE_CHUNKSIZE; \
								if (end_edge_index > (tg)->nglobaledges) end_edge_index = (tg)->nglobaledges; \
								/* fprintf(stderr, "%d trying to read offset = %" PRId64 ", count = %" PRId64 "\n", rank, start_edge_index, edge_count_i); */ \
								MPI_File_read_at_all_end((tg)->edgefile, edge_data_from_file, MPI_STATUS_IGNORE); \
							} else { \
								start_edge_index = int64_min(FILE_CHUNKSIZE * block_idx, (tg)->edgememory_size); \
								end_edge_index = int64_min(start_edge_index + FILE_CHUNKSIZE, (tg)->edgememory_size); \
							} \
							edge_count_i = end_edge_index - start_edge_index; \
							const packed_edge* restrict const user_buf = ((tg)->data_in_file ? edge_data_from_file : (tg)->edgememory + start_edge_index); \
							ptrdiff_t const user_buf_count = edge_count_i; \
							assert (user_buf != NULL); \
							assert (user_buf_count >= 0); \
							assert (tuple_graph_max_bufsize((tg)) >= user_buf_count); \
							int iteration_count = 0; (void)iteration_count; \
							int buffer_released_this_iter = 0; /* To allow explicit buffer release to be optional */ \
							while (1) { \
								/* Prevent continue */ assert (iteration_count == 0); \
								buffer_released_this_iter = 0; \
								{
#endif

#define ITERATE_TUPLE_GRAPH_BLOCK_NUMBER (block_idx)
#define ITERATE_TUPLE_GRAPH_BREAK /* Must be done collectively and before ITERATE_TUPLE_GRAPH_RELEASE_BUFFER */ \
									break_from_block_loop = 1; \
									break
#define ITERATE_TUPLE_GRAPH_RELEASE_BUFFER \
									do { \
										if ((tg)->data_in_file && block_idx + 1 < block_limit) { \
											MPI_Offset start_edge_index = FILE_CHUNKSIZE * (MUL_SIZE((block_idx) + 1) + rank); \
											if (start_edge_index > (tg)->nglobaledges) start_edge_index = (tg)->nglobaledges; \
											edge_count_i = (tg)->nglobaledges - start_edge_index; \
											if (edge_count_i > FILE_CHUNKSIZE) edge_count_i = FILE_CHUNKSIZE; \
											MPI_File_read_at_all_begin((tg)->edgefile, start_edge_index, edge_data_from_file, edge_count_i, packed_edge_mpi_type); \
											buffer_released_this_iter = 1; \
										} \
									} while (0)
#define ITERATE_TUPLE_GRAPH_END \
									if (!buffer_released_this_iter) ITERATE_TUPLE_GRAPH_RELEASE_BUFFER; \
								} \
								if (break_from_block_loop) ITERATE_TUPLE_GRAPH_RELEASE_BUFFER; \
								iteration_count = 1; \
								break; \
							} \
							/* Prevent user break */ assert (iteration_count == 1); \
							if (break_from_block_loop) break; \
						} \
						if (edge_data_from_file) free(edge_data_from_file); \
					} while (0)




					static inline int64_t tuple_graph_max_bufsize(const tuple_graph* tg) {
						return FILE_CHUNKSIZE;
					}

#ifdef __cplusplus
					extern "C" {
#endif

						void setup_globals(void); /* In utils.c */
						void cleanup_globals(void); /* In utils.c */
						int lg_int64_t(int64_t x); /* In utils.c */
						void* xMPI_Alloc_mem(size_t nbytes); /* In utils.c */
						void* xmalloc(size_t nbytes); /* In utils.c */
						void* xcalloc(size_t n, size_t unit); /* In utils.c */

						int validate_result(int isbfs, const tuple_graph* const tg, const size_t nlocalverts, const int64_t root, int64_t* const pred, float * dist, int64_t* const edge_visit_count_ptr); /* In validate.c */

						/* Definitions in each BFS file, using static global variables for internal
						 * storage: */
						void make_graph_data_structure(const tuple_graph* const tg);
						void free_graph_data_structure(void);
						void run_bfs(int64_t root, int64_t* pred);
						void get_edge_count_for_teps(int64_t* edge_visit_count);
						void clean_pred(int64_t* pred);
						size_t get_nlocalverts_for_pred(void);
						/* Definitions in SSSP file in case this kernel is implemented */
#ifdef SSSP
						void run_sssp(int64_t root, int64_t* pred, float * dist_shortest);
						void clean_shortest(float * dist);
#endif

						static inline size_t size_min(size_t a, size_t b) {
							return a < b ? a : b;
						}

						static inline ptrdiff_t ptrdiff_min(ptrdiff_t a, ptrdiff_t b) {
							return a < b ? a : b;
						}

						static inline int64_t int64_min(int64_t a, int64_t b) {
							return a < b ? a : b;
						}

						/* Chunk size for blocks of one-sided operations; a fence is inserted after (at
						 * most) each CHUNKSIZE one-sided operations. */
#define CHUNKSIZE (1 << 22)
#define HALF_CHUNKSIZE ((CHUNKSIZE) / 2)

						/* Bitmap size (in bytes) to keep for each node when determining BFS roots. */
#define BITMAPSIZE (1UL << 29)

#ifdef __cplusplus
					}
#endif

#endif /* COMMON_H */
