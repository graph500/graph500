/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef APPLY_PERMUTATION_MPI_H
#define APPLY_PERMUTATION_MPI_H

#ifdef GRAPH_GENERATOR_MPI

#include <stdint.h>
#include "splittable_mrg.h"
#include "graph_generator.h"
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Internal function to get various information about an uneven block
 * distribution of some data. */
void gather_block_distribution_info(MPI_Comm comm, const int64_t local_block_size, const int64_t global_block_size, int64_t* block_displs /* size = MPI comm size + 1 */, int** block_owner_table_ptr /* malloc'ed in here */, int64_t** block_owner_cutoff_ptr /* malloc'ed in here */, int* lg_minblocksize_ptr, int64_t* maxblocksize_ptr);

/* Internal function to apply a distributed permutation to a distributed set of
 * edges. */
void apply_permutation_mpi(MPI_Comm comm, const int64_t local_perm_size, const int64_t* const local_vertex_perm, const int64_t N, const int64_t nedges, int64_t* result);

#ifdef __cplusplus
}
#endif

#endif /* GRAPH_GENERATOR_MPI */

#endif /* APPLY_PERMUTATION_MPI_H */
