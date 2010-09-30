/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef SCRAMBLE_EDGES_H
#define SCRAMBLE_EDGES_H

#include <stdint.h>
#include "splittable_mrg.h"
#include "graph_generator.h"

#ifdef GRAPH_GENERATOR_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* This version is for sequential machines, OpenMP, and the XMT. */
void scramble_edges_shared(uint64_t userseed1, uint64_t userseed2, int64_t nedges, int64_t* result /* Input and output array of edges (size = 2 * nedges) */);

#ifdef GRAPH_GENERATOR_MPI
/* For MPI distributed memory. */
void scramble_edges_mpi(MPI_Comm comm,
                        const uint64_t userseed1, const uint64_t userseed2,
                        const int64_t local_nedges_in,
                        const int64_t* const local_edges_in,
                        int64_t* const local_nedges_out_ptr,
                        int64_t** const local_edges_out_ptr /* Allocated using xmalloc() by scramble_edges_mpi */);
#endif

#ifdef __cplusplus
}
#endif

#endif /* SCRAMBLE_EDGES_H */
