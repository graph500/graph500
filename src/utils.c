/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <mpi.h>
#include <assert.h>
#include "common.h"

int rank, size;
#ifdef SIZE_MUST_BE_A_POWER_OF_TWO
int lgsize;
#endif
MPI_Datatype packed_edge_mpi_type;

void setup_globals() {
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

#ifdef SIZE_MUST_BE_A_POWER_OF_TWO
	if (/* Check for power of 2 */ (size & (size - 1)) != 0) {
		fprintf(stderr, "Number of processes %d is not a power of two, yet SIZE_MUST_BE_A_POWER_OF_TWO is defined in common.h.\n", size);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	for (lgsize = 0; lgsize < size; ++lgsize) {
		if ((1 << lgsize) == size) break;
	}
	assert (lgsize < size);
#endif

	int blocklengths[] = {1, 1, 1};
	MPI_Aint displs[] = {0, 0, 0};
	packed_edge temp;
	MPI_Aint temp_addr, fld_addr;
	MPI_Get_address(&temp, &temp_addr);
#ifdef GENERATOR_USE_PACKED_EDGE_TYPE
	MPI_Get_address(&temp.v0_low, &fld_addr); displs[0] = fld_addr - temp_addr;
	MPI_Get_address(&temp.v1_low, &fld_addr); displs[1] = fld_addr - temp_addr;
	MPI_Get_address(&temp.high,   &fld_addr); displs[2] = fld_addr - temp_addr;
	MPI_Type_create_hindexed(3, blocklengths, displs, MPI_UINT32_T, &packed_edge_mpi_type);
#else
	MPI_Get_address(&temp.v0, &fld_addr); displs[0] = fld_addr - temp_addr;
	MPI_Get_address(&temp.v1, &fld_addr); displs[1] = fld_addr - temp_addr;
	MPI_Type_create_hindexed(2, blocklengths, displs, MPI_INT64_T, &packed_edge_mpi_type);
#endif
	MPI_Type_commit(&packed_edge_mpi_type);
}

void cleanup_globals(void) {
	MPI_Type_free(&packed_edge_mpi_type);
}

int lg_int64_t(int64_t x) { /* Round up */
	assert (x > 0);
	--x;
	int result = 0;
	while ((x >> result) != 0) ++result;
	return result;
}

void* xMPI_Alloc_mem(size_t nbytes) {
	void* p;
	MPI_Alloc_mem(nbytes, MPI_INFO_NULL, &p);
	if (nbytes != 0 && !p) {
		fprintf(stderr, "MPI_Alloc_mem failed for size %zu\n", nbytes);
		abort();
	}
	return p;
}
