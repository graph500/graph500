/* Copyright (C) 2011 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include "common.h"
#include "onesided.h"
#include <mpi.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

/* One-sided wrapper to allow emulation; a good MPI should be able to handle
 * the version in this file. */

#ifndef EMULATE_ONE_SIDED

/* Gather from one array into another. */
struct gather {
  void* input;
  size_t elt_size;
  void* output;
  MPI_Datatype datatype;
  int valid;
  MPI_Win win;
};

gather* init_gather(void* input, size_t input_count, size_t elt_size, void* output, size_t output_count, size_t nrequests_max, MPI_Datatype dt) {
  gather* g = (gather*)xmalloc(sizeof(gather));
  g->input = input;
  g->elt_size = elt_size;
  g->output = output;
  g->datatype = dt;
  g->valid = 0;
  MPI_Win_create(input, input_count * elt_size, elt_size, MPI_INFO_NULL, MPI_COMM_WORLD, &g->win);
  return g;
}

void destroy_gather(gather* g) {
  assert (!g->valid);
  MPI_Win_free(&g->win);
  free(g);
}

void begin_gather(gather* g) {
  assert (!g->valid);
  g->valid = 1;
  MPI_Win_fence(MPI_MODE_NOPRECEDE | MPI_MODE_NOPUT, g->win);
}

void add_gather_request(gather* g, size_t local_idx, int remote_rank, size_t remote_idx, size_t req_id) {
  assert (g->valid);
#pragma omp critical
  MPI_Get(g->output + local_idx * g->elt_size, 1, g->datatype, remote_rank, remote_idx, 1, g->datatype, g->win);
}

void end_gather(gather* g) {
  assert (g->valid);
  MPI_Win_fence(MPI_MODE_NOSUCCEED, g->win);
  g->valid = 0;
}

/* Scatter a constant to various locations in an array. */
struct scatter_constant {
  void* array;
  size_t elt_size;
  void* constant;
  MPI_Datatype datatype;
  int valid;
  MPI_Win win;
};

scatter_constant* init_scatter_constant(void* array, size_t array_count, size_t elt_size, void* constant, size_t nrequests_max, MPI_Datatype dt) {
  scatter_constant* sc = (scatter_constant*)xmalloc(sizeof(scatter_constant));
  sc->array = array;
  sc->elt_size = elt_size;
  sc->constant = constant;
  sc->datatype = dt;
  sc->valid = 0;
  MPI_Win_create(array, array_count * elt_size, elt_size, MPI_INFO_NULL, MPI_COMM_WORLD, &sc->win);
  return sc;
}

void destroy_scatter_constant(scatter_constant* sc) {
  assert (!sc->valid);
  MPI_Win_free(&sc->win);
  free(sc);
}

void begin_scatter_constant(scatter_constant* sc) {
  assert (!sc->valid);
  sc->valid = 1;
  MPI_Win_fence(MPI_MODE_NOPRECEDE, sc->win);
}

void add_scatter_constant_request(scatter_constant* sc, int remote_rank, size_t remote_idx, size_t req_id) {
  assert (sc->valid);
#pragma omp critical
  MPI_Put(sc->constant, 1, sc->datatype, remote_rank, remote_idx, 1, sc->datatype, sc->win);
}

void end_scatter_constant(scatter_constant* sc) {
  assert (sc->valid);
  MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE, sc->win);
  sc->valid = 0;
}

/* Scatter values to various locations in an array using MPI_REPLACE. */
struct scatter {
  void* array;
  size_t elt_size;
  size_t request_count;
  size_t nrequests_max;
  char* send_data;
  MPI_Datatype datatype;
  int valid;
  MPI_Win win;
};

scatter* init_scatter(void* array, size_t array_count, size_t elt_size, size_t nrequests_max, MPI_Datatype dt) {
  scatter* sc = (scatter*)xmalloc(sizeof(scatter));
  sc->array = array;
  sc->elt_size = elt_size;
  sc->request_count = 0;
  sc->nrequests_max = nrequests_max;
  sc->send_data = xmalloc(nrequests_max * elt_size);
  sc->datatype = dt;
  sc->valid = 0;
  MPI_Win_create(array, array_count * elt_size, elt_size, MPI_INFO_NULL, MPI_COMM_WORLD, &sc->win);
  return sc;
}

void destroy_scatter(scatter* sc) {
  assert (!sc->valid);
  MPI_Win_free(&sc->win);
  free(sc->send_data);
  free(sc);
}

void begin_scatter(scatter* sc) {
  assert (!sc->valid);
  sc->valid = 1;
  sc->request_count = 0;
  MPI_Win_fence(MPI_MODE_NOPRECEDE, sc->win);
}

void add_scatter_request(scatter* sc, const char* local_data, int remote_rank, size_t remote_idx, size_t req_id) {
  assert (sc->valid);
  assert (sc->request_count < sc->nrequests_max);
  memcpy(sc->send_data + sc->request_count * sc->elt_size, local_data, sc->elt_size);
#pragma omp critical
  MPI_Put(sc->send_data + sc->request_count * sc->elt_size, 1, sc->datatype, remote_rank, remote_idx, 1, sc->datatype, sc->win);
  ++sc->request_count;
}

void end_scatter(scatter* sc) {
  assert (sc->valid);
  MPI_Win_fence(MPI_MODE_NOSUCCEED | MPI_MODE_NOSTORE, sc->win);
  sc->valid = 0;
}

#endif /* !EMULATE_ONE_SIDED */
