/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef PERMUTATION_GEN_H
#define PERMUTATION_GEN_H

#include <stdint.h>
#include "splittable_mrg.h"
#include "graph_generator.h"

#ifdef GRAPH_GENERATOR_MPI
#include <mpi.h>
#endif

/* This version is for sequential machines and the XMT. */
void rand_sort_shared(mrg_state* st, int64_t n, int64_t* result /* Array of size n */);

#ifdef GRAPH_GENERATOR_MPI
/* For MPI distributed memory. */
void rand_sort_mpi(MPI_Comm comm, mrg_state* st, int64_t n,
                   int64_t* result_size_ptr,
                   int64_t** result_ptr /* Allocated using malloc() by
                   rand_sort_mpi(), must be free()d by user */);
#endif

#endif /* PERMUTATION_GEN_H */
