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

/* One-sided emulation since many MPI implementations don't have good
 * performance and/or fail when using many one-sided operations between fences.
 * Only the necessary operations are implemented, and only MPI_MODE_NOPRECEDE
 * and MPI_MODE_NOSUCCEED fences can be used. */

#ifdef EMULATE_ONE_SIDED

/* Gather from one array into another. */
struct gather {
  void* input;
  size_t input_count;
  size_t elt_size;
  void* output;
  size_t output_count;
  size_t nrequests_max;
  MPI_Datatype datatype;
  int valid;
  MPI_Comm comm;
  size_t* local_indices;
  int* remote_ranks;
  MPI_Aint* remote_indices;
  int comm_size;
  int* send_counts;
  int* send_offsets;
  int* recv_counts;
  int* recv_offsets;
};

gather* init_gather(void* input, size_t input_count, size_t elt_size, void* output, size_t output_count, size_t nrequests_max, MPI_Datatype dt) {
  gather* g = (gather*)xmalloc(sizeof(gather));
  g->input = input;
  g->input_count = input_count;
  g->elt_size = elt_size;
  g->output = output;
  g->output_count = output_count;
  g->nrequests_max = nrequests_max;
  g->datatype = dt;
  g->valid = 0;
  MPI_Comm_dup(MPI_COMM_WORLD, &g->comm);
  g->local_indices = (size_t*)xmalloc(nrequests_max * sizeof(size_t));
  g->remote_ranks = (int*)xmalloc(nrequests_max * sizeof(int));
  g->remote_indices = (MPI_Aint*)xmalloc(nrequests_max * sizeof(MPI_Aint));
  MPI_Comm_size(g->comm, &g->comm_size);
  int size = g->comm_size;
  g->send_counts = (int*)xmalloc((size + 1) * sizeof(int));
  g->send_offsets = (int*)xmalloc((size + 2) * sizeof(int));
  g->recv_counts = (int*)xmalloc(size * sizeof(int));
  g->recv_offsets = (int*)xmalloc((size + 1) * sizeof(int));
#ifndef NDEBUG
  int datasize;
  MPI_Type_size(dt, &datasize);
  assert (datasize == (int)elt_size);
#endif
  return g;
}

void destroy_gather(gather* g) {
  assert (!g->valid);
  free(g->local_indices); g->local_indices = NULL;
  free(g->remote_ranks); g->remote_ranks = NULL;
  free(g->remote_indices); g->remote_indices = NULL;
  MPI_Comm_free(&g->comm);
  free(g->send_counts); g->send_counts = NULL;
  free(g->send_offsets); g->send_offsets = NULL;
  free(g->recv_counts); g->recv_counts = NULL;
  free(g->recv_offsets); g->recv_offsets = NULL;
  g->input = NULL;
  g->output = NULL;
  g->input_count = 0;
  g->output_count = 0;
  free(g);
}

void begin_gather(gather* g) {
  assert (!g->valid);
  {size_t i, nr = g->nrequests_max; for (i = 0; i < nr; ++i) g->remote_ranks[i] = g->comm_size;}
  g->valid = 1;
}

void add_gather_request(gather* g, size_t local_idx, int remote_rank, size_t remote_idx, size_t req_id) {
  assert (g->valid);
  assert (remote_rank >= 0 && remote_rank < g->comm_size);
  assert (req_id < g->nrequests_max);
  assert (local_idx < g->output_count);
  g->local_indices[req_id] = local_idx;
  assert (g->remote_ranks[req_id] == g->comm_size);
  g->remote_ranks[req_id] = remote_rank;
  g->remote_indices[req_id] = (MPI_Aint)remote_idx;
}

/* Adapted from histogram_sort_inplace in boost/graph/detail/histogram_sort.hpp
 * */
void histogram_sort_size_tMPI_Aint
       (int* restrict keys,
        const int* restrict rowstart,
        int numkeys,
        size_t* restrict values1,
        MPI_Aint* restrict values2) {
  int* restrict insert_positions = (int*)xmalloc(numkeys * sizeof(int));
  memcpy(insert_positions, rowstart, numkeys * sizeof(int));
  int i;
  for (i = 0; i < rowstart[numkeys]; ++i) {
    int key = keys[i];
    assert (key >= 0 && key < numkeys);
    // fprintf(stderr, "i = %d, key = %d\n", i, key);
    while (!(i >= rowstart[key] && i < insert_positions[key])) {
      int target_pos = insert_positions[key]++;
      // fprintf(stderr, "target_pos = %d from limit %d\n", target_pos, rowstart[key + 1]);
      if (target_pos == i) continue;
      assert (target_pos < rowstart[key + 1]);
      // fprintf(stderr, "swapping [%d] = (%d, %d, %d) with [%d] = (%d, %d, %d)\n", i, keys[i], (int)values1[i], (int)values2[i], target_pos, keys[target_pos], (int)values1[target_pos], (int)values2[target_pos]);
      {int t = keys[i]; key = keys[target_pos]; keys[i] = key; keys[target_pos] = t;}
      {size_t t = values1[i]; values1[i] = values1[target_pos]; values1[target_pos] = t;}
      {MPI_Aint t = values2[i]; values2[i] = values2[target_pos]; values2[target_pos] = t;}
      assert (key >= 0 && key < numkeys);
    }
    // fprintf(stderr, "done\n");
  }
  for (i = 1; i < rowstart[numkeys]; ++i) {
    assert (keys[i] >= keys[i - 1]);
  }
  free(insert_positions);
}

void end_gather(gather* g) {
  assert (g->valid);
  int size = g->comm_size;
  int* restrict send_counts = g->send_counts;
  int* restrict send_offsets = g->send_offsets;
  int* restrict recv_counts = g->recv_counts;
  int* restrict recv_offsets = g->recv_offsets;
  size_t* restrict local_indices = g->local_indices;
  int* restrict remote_ranks = g->remote_ranks;
  MPI_Aint* restrict remote_indices = g->remote_indices;
  const char* restrict input = (const char*)g->input;
  char* restrict output = (char*)g->output;
  size_t elt_size = g->elt_size;
  MPI_Comm comm = g->comm;
  MPI_Datatype datatype = g->datatype;
  size_t nrequests_max = g->nrequests_max;
#ifndef NDEBUG
  size_t input_count = g->input_count;
#endif
  memset(send_counts, 0, (size + 1) * sizeof(int));
  size_t i;
  for (i = 0; i < nrequests_max; ++i) {
    assert (remote_ranks[i] >= 0 && remote_ranks[i] < size + 1);
    ++send_counts[remote_ranks[i]];
  }
  send_offsets[0] = 0;
  for (i = 0; i < (size_t)size + 1; ++i) {
    assert (send_counts[i] >= 0);
    send_offsets[i + 1] = send_offsets[i] + send_counts[i];
  }
  assert (send_offsets[size + 1] == (int)nrequests_max);
  histogram_sort_size_tMPI_Aint(remote_ranks, send_offsets, size + 1, local_indices, remote_indices);
  assert (send_offsets[size] == send_offsets[size + 1] || remote_ranks[send_offsets[size]] == size);
  MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);
  recv_offsets[0] = 0;
  for (i = 0; i < (size_t)size; ++i) {
    assert (recv_counts[i] >= 0);
    recv_offsets[i + 1] = recv_offsets[i] + recv_counts[i];
  }
  MPI_Aint* restrict recv_data = (MPI_Aint*)xmalloc(recv_offsets[size] * sizeof(MPI_Aint));
  MPI_Alltoallv(remote_indices, send_counts, send_offsets, MPI_AINT, recv_data, recv_counts, recv_offsets, MPI_AINT, comm);
  char* restrict reply_data = xmalloc(recv_offsets[size] * elt_size);
  for (i = 0; i < (size_t)recv_offsets[size]; ++i) {
    assert (recv_data[i] >= 0 && recv_data[i] < input_count);
    memcpy(reply_data + i * elt_size, input + recv_data[i] * elt_size, elt_size);
  }
  free(recv_data);
  char* restrict recv_reply_data = xmalloc(send_offsets[size] * elt_size);
  MPI_Alltoallv(reply_data, recv_counts, recv_offsets, datatype, recv_reply_data, send_counts, send_offsets, datatype, comm);
  free(reply_data);
  for (i = 0; i < nrequests_max; ++i) {
    if (remote_ranks[i] >= 0 && remote_ranks[i] < size) {
      memcpy(output + local_indices[i] * elt_size, recv_reply_data + i * elt_size, elt_size);
    }
  }
  free(recv_reply_data);
  g->valid = 0;
}

struct scatter_constant {
  void* array;
  size_t array_count;
  size_t elt_size;
  void* constant;
  size_t nrequests_max;
  int valid;
  MPI_Comm comm;
  int* remote_ranks;
  MPI_Aint* remote_indices;
  int comm_size;
  int* send_counts;
  int* send_offsets;
  int* recv_counts;
  int* recv_offsets;
};

scatter_constant* init_scatter_constant(void* array, size_t array_count, size_t elt_size, void* constant, size_t nrequests_max, MPI_Datatype dt /* unused */) {
  (void)dt;
  scatter_constant* sc = (scatter_constant*)xmalloc(sizeof(scatter_constant));
  sc->array = array;
  sc->array_count = array_count;
  sc->elt_size = elt_size;
  sc->constant = constant;
  sc->nrequests_max = nrequests_max;
  sc->valid = 0;
  MPI_Comm_dup(MPI_COMM_WORLD, &sc->comm);
  sc->remote_ranks = (int*)xmalloc(nrequests_max * sizeof(int));
  sc->remote_indices = (MPI_Aint*)xmalloc(nrequests_max * sizeof(MPI_Aint));
  MPI_Comm_size(sc->comm, &sc->comm_size);
  int size = sc->comm_size;
  sc->send_counts = (int*)xmalloc((size + 1) * sizeof(int));
  sc->send_offsets = (int*)xmalloc((size + 2) * sizeof(int));
  sc->recv_counts = (int*)xmalloc(size * sizeof(int));
  sc->recv_offsets = (int*)xmalloc((size + 1) * sizeof(int));
  return sc;
}

void destroy_scatter_constant(scatter_constant* sc) {
  assert (!sc->valid);
  free(sc->remote_ranks); sc->remote_ranks = NULL;
  free(sc->remote_indices); sc->remote_indices = NULL;
  MPI_Comm_free(&sc->comm);
  free(sc->send_counts); sc->send_counts = NULL;
  free(sc->send_offsets); sc->send_offsets = NULL;
  free(sc->recv_counts); sc->recv_counts = NULL;
  free(sc->recv_offsets); sc->recv_offsets = NULL;
  free(sc);
}

void begin_scatter_constant(scatter_constant* sc) {
  assert (!sc->valid);
  {size_t i, nr = sc->nrequests_max; for (i = 0; i < nr; ++i) sc->remote_ranks[i] = sc->comm_size;}
  sc->valid = 1;
}

void add_scatter_constant_request(scatter_constant* sc, int remote_rank, size_t remote_idx, size_t req_id) {
  assert (sc->valid);
  assert (remote_rank >= 0 && remote_rank < sc->comm_size);
  assert (req_id < sc->nrequests_max);
  assert (sc->remote_ranks[req_id] == sc->comm_size);
  sc->remote_ranks[req_id] = remote_rank;
  sc->remote_indices[req_id] = (MPI_Aint)remote_idx;
}

/* Adapted from histogram_sort_inplace in boost/graph/detail/histogram_sort.hpp
 * */
void histogram_sort_MPI_Aint
       (int* restrict keys,
        const int* restrict rowstart,
        int numkeys,
        MPI_Aint* restrict values1) {
  int* restrict insert_positions = (int*)xmalloc(numkeys * sizeof(int));
  memcpy(insert_positions, rowstart, numkeys * sizeof(int));
  int i;
  for (i = 0; i < rowstart[numkeys]; ++i) {
    int key = keys[i];
    assert (key >= 0 && key < numkeys);
    while (!(i >= rowstart[key] && i < insert_positions[key])) {
      int target_pos = insert_positions[key]++;
      if (target_pos == i) continue;
      assert (target_pos < rowstart[key + 1]);
      {int t = keys[i]; key = keys[target_pos]; keys[i] = key; keys[target_pos] = t;}
      {MPI_Aint t = values1[i]; values1[i] = values1[target_pos]; values1[target_pos] = t;}
      assert (key >= 0 && key < numkeys);
    }
  }
  for (i = 1; i < rowstart[numkeys]; ++i) {
    assert (keys[i] >= keys[i - 1]);
  }
  free(insert_positions);
}

void end_scatter_constant(scatter_constant* sc) {
  assert (sc->valid);
  int size = sc->comm_size;
  int* restrict send_counts = sc->send_counts;
  int* restrict send_offsets = sc->send_offsets;
  int* restrict recv_counts = sc->recv_counts;
  int* restrict recv_offsets = sc->recv_offsets;
  int* restrict remote_ranks = sc->remote_ranks;
  MPI_Aint* restrict remote_indices = sc->remote_indices;
  char* restrict array = (char*)sc->array;
  const char* restrict constant = (const char*)sc->constant;
  size_t elt_size = sc->elt_size;
  MPI_Comm comm = sc->comm;
  size_t nrequests_max = sc->nrequests_max;
#ifndef NDEBUG
  size_t array_count = sc->array_count;
#endif
  memset(send_counts, 0, (size + 1) * sizeof(int));
  size_t i;
  for (i = 0; i < nrequests_max; ++i) {
    assert (remote_ranks[i] >= 0 && remote_ranks[i] < size + 1);
    ++send_counts[remote_ranks[i]];
  }
  send_offsets[0] = 0;
  for (i = 0; i < (size_t)size + 1; ++i) {
    send_offsets[i + 1] = send_offsets[i] + send_counts[i];
  }
  histogram_sort_MPI_Aint(remote_ranks, send_offsets, size + 1, remote_indices);
  MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);
  recv_offsets[0] = 0;
  for (i = 0; i < (size_t)size; ++i) {
    recv_offsets[i + 1] = recv_offsets[i] + recv_counts[i];
  }
  MPI_Aint* restrict recv_data = (MPI_Aint*)xmalloc(recv_offsets[size] * sizeof(MPI_Aint));
  MPI_Alltoallv(remote_indices, send_counts, send_offsets, MPI_AINT, recv_data, recv_counts, recv_offsets, MPI_AINT, comm);
  for (i = 0; i < (size_t)recv_offsets[size]; ++i) {
    assert (recv_data[i] >= 0 && recv_data[i] < array_count);
    memcpy(array + recv_data[i] * elt_size, constant, elt_size);
  }
  free(recv_data);
  sc->valid = 0;
}

struct scatter {
  void* array;
  size_t array_count;
  size_t elt_size;
  char* send_data;
  size_t nrequests_max;
  MPI_Datatype datatype;
  int valid;
  MPI_Comm comm;
  int* remote_ranks;
  MPI_Aint* remote_indices;
  int comm_size;
  int* send_counts;
  int* send_offsets;
  int* recv_counts;
  int* recv_offsets;
};

scatter* init_scatter(void* array, size_t array_count, size_t elt_size, size_t nrequests_max, MPI_Datatype dt) {
  scatter* sc = (scatter*)xmalloc(sizeof(scatter));
  sc->array = array;
  sc->array_count = array_count;
  sc->elt_size = elt_size;
  sc->send_data = xmalloc(nrequests_max * elt_size);
  sc->nrequests_max = nrequests_max;
  sc->datatype = dt;
  sc->valid = 0;
  MPI_Comm_dup(MPI_COMM_WORLD, &sc->comm);
  sc->remote_ranks = (int*)xmalloc(nrequests_max * sizeof(int));
  sc->remote_indices = (MPI_Aint*)xmalloc(nrequests_max * sizeof(MPI_Aint));
  MPI_Comm_size(sc->comm, &sc->comm_size);
  int size = sc->comm_size;
  sc->send_counts = (int*)xmalloc((size + 1) * sizeof(int));
  sc->send_offsets = (int*)xmalloc((size + 2) * sizeof(int));
  sc->recv_counts = (int*)xmalloc(size * sizeof(int));
  sc->recv_offsets = (int*)xmalloc((size + 1) * sizeof(int));
  return sc;
}

void destroy_scatter(scatter* sc) {
  assert (!sc->valid);
  free(sc->send_data); sc->send_data = NULL;
  free(sc->remote_ranks); sc->remote_ranks = NULL;
  free(sc->remote_indices); sc->remote_indices = NULL;
  MPI_Comm_free(&sc->comm);
  free(sc->send_counts); sc->send_counts = NULL;
  free(sc->send_offsets); sc->send_offsets = NULL;
  free(sc->recv_counts); sc->recv_counts = NULL;
  free(sc->recv_offsets); sc->recv_offsets = NULL;
  free(sc);
}

void begin_scatter(scatter* sc) {
  assert (!sc->valid);
  {size_t i, nr = sc->nrequests_max; for (i = 0; i < nr; ++i) sc->remote_ranks[i] = sc->comm_size;}
  sc->valid = 1;
}

void add_scatter_request(scatter* sc, const char* local_data, int remote_rank, size_t remote_idx, size_t req_id) {
  assert (sc->valid);
  assert (req_id < sc->nrequests_max);
  assert (remote_rank >= 0 && remote_rank < sc->comm_size);
  memcpy(sc->send_data + req_id * sc->elt_size, local_data, sc->elt_size);
  assert (sc->remote_ranks[req_id] == sc->comm_size);
  sc->remote_ranks[req_id] = remote_rank;
  sc->remote_indices[req_id] = (MPI_Aint)remote_idx;
}

/* Adapted from histogram_sort_inplace in boost/graph/detail/histogram_sort.hpp
 * */
void histogram_sort_MPI_Aintcharblock
       (int* restrict keys,
        const int* restrict rowstart,
        int numkeys,
        MPI_Aint* restrict values1,
        char* restrict values2,
        size_t elt_size2) {
  int* restrict insert_positions = (int*)xmalloc(numkeys * sizeof(int));
  memcpy(insert_positions, rowstart, numkeys * sizeof(int));
  int i;
  for (i = 0; i < rowstart[numkeys]; ++i) {
    int key = keys[i];
    assert (key >= 0 && key < numkeys);
    fprintf(stderr, "i = %d, key = %d\n", i, key);
    while (!(i >= rowstart[key] && i < insert_positions[key])) {
      int target_pos = insert_positions[key]++;
      fprintf(stderr, "target_pos = %d from limit %d\n", target_pos, rowstart[key + 1]);
      if (target_pos == i) continue;
      assert (target_pos < rowstart[key + 1]);
      {int t = keys[i]; key = keys[target_pos]; keys[i] = key; keys[target_pos] = t;}
      {MPI_Aint t = values1[i]; values1[i] = values1[target_pos]; values1[target_pos] = t;}
      {char t[elt_size2]; memcpy(t, values2 + i * elt_size2, elt_size2); memcpy(values2 + i * elt_size2, values2 + target_pos * elt_size2, elt_size2); memcpy(values2 + target_pos * elt_size2, t, elt_size2);}
      assert (key >= 0 && key < numkeys);
    }
    fprintf(stderr, "done\n");
  }
  for (i = 1; i < rowstart[numkeys]; ++i) {
    assert (keys[i] >= keys[i - 1]);
  }
  free(insert_positions);
}

void end_scatter(scatter* sc) {
  assert (sc->valid);
  int size = sc->comm_size;
  int* restrict send_counts = sc->send_counts;
  int* restrict send_offsets = sc->send_offsets;
  int* restrict recv_counts = sc->recv_counts;
  int* restrict recv_offsets = sc->recv_offsets;
  char* restrict send_data = sc->send_data;
  int* restrict remote_ranks = sc->remote_ranks;
  MPI_Aint* restrict remote_indices = sc->remote_indices;
  char* restrict array = (char*)sc->array;
  size_t elt_size = sc->elt_size;
  MPI_Comm comm = sc->comm;
  size_t nrequests_max = sc->nrequests_max;
#ifndef NDEBUG
  size_t array_count = sc->array_count;
#endif
  memset(send_counts, 0, (size + 1) * sizeof(int));
  size_t i;
  for (i = 0; i < nrequests_max; ++i) {
    assert (remote_ranks[i] >= 0 && remote_ranks[i] < size + 1);
    ++send_counts[remote_ranks[i]];
  }
  send_offsets[0] = 0;
  for (i = 0; i < (size_t)size + 1; ++i) {
    send_offsets[i + 1] = send_offsets[i] + send_counts[i];
  }
  histogram_sort_MPI_Aintcharblock(remote_ranks, send_offsets, size + 1, remote_indices, send_data, elt_size);
  MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);
  recv_offsets[0] = 0;
  for (i = 0; i < (size_t)size; ++i) {
    recv_offsets[i + 1] = recv_offsets[i] + recv_counts[i];
  }
  MPI_Aint* restrict recv_indices = (MPI_Aint*)xmalloc(recv_offsets[size] * sizeof(MPI_Aint));
  char* restrict recv_data = (char*)xmalloc(recv_offsets[size] * elt_size);
  MPI_Alltoallv(remote_indices, send_counts, send_offsets, MPI_AINT, recv_indices, recv_counts, recv_offsets, MPI_AINT, comm);
  MPI_Alltoallv(send_data, send_counts, send_offsets, sc->datatype, recv_data, recv_counts, recv_offsets, sc->datatype, comm);
  for (i = 0; i < (size_t)recv_offsets[size]; ++i) {
    assert (recv_indices[i] >= 0 && recv_indices[i] < array_count);
    memcpy(array + recv_indices[i] * elt_size, recv_data + i * elt_size, elt_size);
  }
  free(recv_data);
  free(recv_indices);
  sc->valid = 0;
}

#endif /* EMULATE_ONE_SIDED */
