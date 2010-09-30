/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef GRAPH_GENERATOR_H
#define GRAPH_GENERATOR_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Settings for user modification ----------------------------------- */

#ifndef MODIFY_PARAMS_AT_EACH_LEVEL
/* Add noise into each subblock's Kronecker parameters (as is done in SSCA #2,
 * but here the same mutation is used for all edges within a block) */
/* #define MODIFY_PARAMS_AT_EACH_LEVEL */
#endif

#define INT64_T_MPI_TYPE MPI_LONG_LONG /* Should be MPI_INT64_T */

/* End of user settings ----------------------------------- */

/* Graph generator #define settings -- set with benchmark configuration, but
 * should not be changed by individual runners of the benchmark (other than
 * GRAPHGEN_DISTRIBUTED_MEMORY). */

#ifndef GRAPHGEN_SETTINGS_DEFINED
#define GRAPHGEN_SETTINGS_DEFINED

/* Define if output array is local to each rank (and the output array pointer
 * points to the beginning of the local part) instead of being global (i.e.,
 * having the output array pointer pointing to the same place on all ranks); in
 * either case, the output array is still broken into independent pieces). */
/* #define GRAPHGEN_DISTRIBUTED_MEMORY -- Set this in your Makefile */

/* Define if Kronecker parameters should be modified within each sub-block of
 * the adjacency matrix (like SSCA #2 does). */
/* #define GRAPHGEN_MODIFY_PARAMS_AT_EACH_LEVEL */

/* Size of Kronecker initiator matrix. */
#define GRAPHGEN_INITIATOR_SIZE 2

/* Define to clip edges to one half of adjacency matrix to create undirected
 * graphs. */
#define GRAPHGEN_UNDIRECTED

/* Define to return graph edges in a struct that contains multiplicities,
 * rather than just as an array of endpoints with duplicates and self-loops
 * removed. */
/* #define GRAPHGEN_KEEP_MULTIPLICITIES */

/* Define to keep self-loops in output array rather than marking them as unused
 * slots. */
#define GRAPHGEN_KEEP_SELF_LOOPS

/* Define to keep duplicate edges in output array rather than marking them as
 * unused slots. */
#define GRAPHGEN_KEEP_DUPLICATES

#endif /* GRAPHGEN_SETTINGS_DEFINED */

#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
typedef struct generated_edge {
  int64_t src;
  int64_t tgt;
  int64_t multiplicity;
} generated_edge;
#endif

int64_t compute_edge_array_size(
       int rank, int size,
       int64_t M);

void generate_kronecker(
       int rank, int size,
       const uint_fast32_t seed[5] /* All values in [0, 2^31 - 1) */,
       int logN /* In base initiator_size */,
       int64_t M,
       const double initiator[ /* initiator_size * initiator_size */ ],
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
       generated_edge* const edges /* Size >= compute_edge_array_size(rank,
       size, M), must be zero-initialized; points to beginning of local chunk
       when output_array_is_local is 1 and beginning of global array when
       output_array_is_local is 0 */
#else
       int64_t* const edges /* Size >= 2 * compute_edge_array_size(rank, size,
       M); two endpoints per edge (= -1 when slot is unused); points to
       beginning of local chunk when output_array_is_local is 1 and beginning
       of global array when output_array_is_local is 0 */
#endif
);

#ifdef __cplusplus
}
#endif

#endif /* GRAPH_GENERATOR_H */
