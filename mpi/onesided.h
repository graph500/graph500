/* Copyright (C) 2011 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include "common.h"
#include <mpi.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#ifndef ONESIDED_H
#define ONESIDED_H

/* Gather from one array into another. */
typedef struct gather gather;

gather* init_gather(void* input, size_t input_count, size_t elt_size, void* output, size_t output_count, size_t nrequests_max, MPI_Datatype dt);
void destroy_gather(gather* g);
void begin_gather(gather* g);
void add_gather_request(gather* g, size_t local_idx, int remote_rank, size_t remote_idx, size_t req_id);
void end_gather(gather* g);

/* Scatter a constant to various locations in an array. */
typedef struct scatter_constant scatter_constant;

scatter_constant* init_scatter_constant(void* array, size_t array_size, size_t elt_size, void* constant, size_t nrequests_max, MPI_Datatype dt);
void destroy_scatter_constant(scatter_constant* sc);
void begin_scatter_constant(scatter_constant* sc);
void add_scatter_constant_request(scatter_constant* sc, int remote_rank, size_t remote_idx, size_t req_id);
void end_scatter_constant(scatter_constant* sc);

/* Scatter values to various locations in an array (this must use MPI_REPLACE).
 * */
typedef struct scatter scatter;

scatter* init_scatter(void* array, size_t array_size, size_t elt_size, size_t nrequests_max, MPI_Datatype dt);
void destroy_scatter(scatter* sc);
void begin_scatter(scatter* sc);
void add_scatter_request(scatter* sc, const char* local_data, int remote_rank, size_t remote_idx, size_t req_id);
void end_scatter(scatter* sc);

/* There are two implementations of these functions, in onesided.c and
 * onesided_emul.c.  The EMULATE_ONE_SIDED #define controls which of those
 * files is used. */

#endif /* ONESIDED_H */

