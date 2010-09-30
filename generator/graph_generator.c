/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "splittable_mrg.h"
#include "btrd_binomial_distribution.h"
#include "graph_generator.h"

#define GRAPHGEN_INITIATOR_SIZE2 (GRAPHGEN_INITIATOR_SIZE * GRAPHGEN_INITIATOR_SIZE)

typedef struct generator_settings { /* Internal settings */
  double initiator[GRAPHGEN_INITIATOR_SIZE2];
  int64_t my_first_edge, my_last_edge;
  int64_t total_nverts;
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  generated_edge* out; /* Start of local (when GRAPHGEN_DISTRIBUTED_MEMORY is defined) or global (when GRAPHGEN_DISTRIBUTED_MEMORY is undefined) output edge list */
#else
  int64_t* out; /* Start of local (when GRAPHGEN_DISTRIBUTED_MEMORY is defined) or global (when GRAPHGEN_DISTRIBUTED_MEMORY is undefined) output edge list */
#endif
} generator_settings;

static int generate_nway_bernoulli(const generator_settings* s, mrg_state* st) {
  double random_number;
  int j;
  random_number = mrg_get_double_orig(st);
  for (j = 0; j + 1 < GRAPHGEN_INITIATOR_SIZE2; ++j) {
    double ini = s->initiator[j];
    if (random_number < ini) {
      return j;
    }
    random_number -= ini;
  }
  return GRAPHGEN_INITIATOR_SIZE2 - 1;
}

static
void make_square_counts(int64_t num_edges, mrg_state* st, const generator_settings* s, int64_t* square_counts) {
  /* fprintf(stderr, "make_square_counts %zu with (%lf, %lf, %lf, %lf)\n", num_edges, a, b, c, d); */
  if (num_edges <= 20) {
    int i;
    int64_t ii;
    for (i = 0; i < GRAPHGEN_INITIATOR_SIZE2; ++i) {
      square_counts[i] = 0;
    }
    for (ii = 0; ii < num_edges; ++ii) {
      /* ++square_counts[generate_nway_bernoulli(s, st)]; -- replaced by code below */
      int j;
      double random_number = mrg_get_double_orig(st);
      for (j = 0; j < GRAPHGEN_INITIATOR_SIZE2; ++j) {
        double ini = s->initiator[j];
        if (random_number < ini || j == GRAPHGEN_INITIATOR_SIZE2 - 1) {
          ++square_counts[j];
          break;
        }
        random_number -= ini;
      }
    }
  } else {
    int64_t num_edges_left = num_edges;
    double divisor = 1.;
    int i;
    for (i = 0; i < GRAPHGEN_INITIATOR_SIZE2 - 1; ++i) {
      square_counts[i] = btrd_binomial_distribution(num_edges_left, s->initiator[i] / divisor, st);
      num_edges_left -= square_counts[i];
      divisor -= s->initiator[i];
    }
    square_counts[GRAPHGEN_INITIATOR_SIZE2 - 1] = num_edges_left;
  }
  /* fprintf(stderr, "--> (%zu, %zu, %zu, %zu)\n", square_counts[0], square_counts[1], square_counts[2], square_counts[3]); */
}

#ifdef GRAPHGEN_MODIFY_PARAMS_AT_EACH_LEVEL
static void alter_params(generator_settings* s, const mrg_transition_matrix* trans, mrg_state* st) {
  double divisor = 0.;
  for (int i = 0; i < GRAPHGEN_INITIATOR_SIZE2; ++i) {
    s->initiator[i] *= .95 + .1 * mrg_get_double(trans, st);
    divisor += s->initiator[i];
  }
  for (int i = 0; i < GRAPHGEN_INITIATOR_SIZE2; ++i) {
    s->initiator[i] /= divisor;
  }
}
#endif

static
void make_one_edge(int64_t base_src, int64_t base_tgt, int64_t nverts, mrg_state* st, const generator_settings* s,
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
                   generated_edge* result
#else
                   int64_t* result
#endif
                   ) {
#ifdef GRAPHGEN_MODIFY_PARAMS_AT_EACH_LEVEL
  generator_settings my_settings_data = s;
  const generator_settings* my_settings = &my_settings_data;
#else
  const generator_settings* my_settings = s;
#endif
  while (nverts > 1) {
    int square = generate_nway_bernoulli(my_settings, st);
    int src_offset = square / GRAPHGEN_INITIATOR_SIZE;
    int tgt_offset = square % GRAPHGEN_INITIATOR_SIZE;
#ifdef GRAPHGEN_UNDIRECTED
    assert (base_src <= base_tgt);
    if (base_src == base_tgt) {
      /* Clip-and-flip for undirected graph */
      if (src_offset > tgt_offset) {
        int temp = src_offset;
        src_offset = tgt_offset;
        tgt_offset = temp;
      }
    }
#endif /* GRAPHGEN_UNDIRECTED */
    nverts /= GRAPHGEN_INITIATOR_SIZE;
    base_src += nverts * src_offset;
    base_tgt += nverts * tgt_offset;
#ifdef GRAPHGEN_MODIFY_PARAMS_AT_EACH_LEVEL
    alter_params(my_settings, trans, st);
#endif
  }
#ifdef GRAPHGEN_KEEP_SELF_LOOPS
  int is_self_loop_to_skip = 0;
#else
  int is_self_loop_to_skip = (base_src == base_tgt);
#endif
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
  result->src = base_src;
  result->tgt = base_tgt;
  result->multiplicity = is_self_loop_to_skip ? 0 : 1;
#else
  result[0] = is_self_loop_to_skip ? -1 : base_src;
  result[1] = is_self_loop_to_skip ? -1 : base_tgt;
#endif
}

static void generate_kronecker_internal(
  const mrg_state* orig_state,
  int64_t first_edge_index,
  int64_t num_edges,
  int64_t nverts,
  const generator_settings* s,
  int64_t base_src,
  int64_t base_tgt) {
  mrg_state state = *orig_state;
  mrg_skip(&state, 0, (base_src + s->total_nverts) / nverts, (base_tgt + s->total_nverts) / nverts);
  int64_t my_first_edge = s->my_first_edge;
  int64_t my_last_edge = s->my_last_edge;
#ifdef GRAPHGEN_UNDIRECTED
  assert (base_src <= base_tgt);
#endif /* GRAPHGEN_UNDIRECTED */
  if (nverts == 1) {
    assert (num_edges != 0);
#ifdef GRAPHGEN_KEEP_SELF_LOOPS
    int is_self_loop_to_skip = 0;
#else
    int is_self_loop_to_skip = (base_src == base_tgt);
#endif
    int i;
    for (i = 0; i < num_edges; ++i) {
      /* Write all edges, filling all slots except the first with edges marked
       * as removed duplicates; the complexity of the loop here is to deal with
       * cases of overflows between nodes (i.e., some of the copies of an edge
       * on one node and some on another).  Also, if the edge is a self-loop
       * and those are being removed, write -1 (or a multiplicity of 0) instead
       * so it will be removed. */
      if (first_edge_index + i >= my_first_edge && first_edge_index + i < my_last_edge) {
#ifdef GRAPHGEN_DISTRIBUTED_MEMORY
        int64_t write_offset = first_edge_index + i - my_first_edge;
#else
        int64_t write_offset = first_edge_index + i;
#endif
#ifdef GRAPHGEN_KEEP_DUPLICATES
        int is_duplicate_to_skip = 0;
#else
        int is_duplicate_to_skip = (i != 0);
#endif
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
        {
#ifdef GRAPHGEN_KEEP_DUPLICATES
          int64_t multiplicity_to_write = 1; /* Write all edges as 1 */
#else
          int64_t multiplicity_to_write = num_edges; /* Write first edge as num_edges and rest as 0 */
#endif
          generated_edge* out_loc = &(s->out[write_offset]);
          out_loc->src = base_src;
          out_loc->tgt = base_tgt;
          out_loc->multiplicity = ((!is_duplicate_to_skip && !is_self_loop_to_skip) ? multiplicity_to_write : 0);
        }
#else
        {
          int64_t* out_loc = &(s->out[2 * write_offset]);
          out_loc[0] = ((!is_duplicate_to_skip && !is_self_loop_to_skip) ? base_src : -1);
          out_loc[1] = ((!is_duplicate_to_skip && !is_self_loop_to_skip) ? base_tgt : -1);
        }
#endif
      }
    }
  } else if (num_edges == 1) {
    if (first_edge_index >= my_first_edge && first_edge_index < my_last_edge) {
#ifdef GRAPHGEN_DISTRIBUTED_MEMORY
      int64_t write_offset = first_edge_index - my_first_edge;
#else
      int64_t write_offset = first_edge_index;
#endif
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
      generated_edge* out_loc = &(s->out[write_offset]);
#else
      int64_t* out_loc = &(s->out[2 * write_offset]);
#endif
      make_one_edge(base_src, base_tgt, nverts, &state, s, out_loc);
    }
  } else {
    int64_t new_nverts;
#ifdef GRAPHGEN_MODIFY_PARAMS_AT_EACH_LEVEL
    generator_settings new_settings_data;
    const generator_settings* new_settings = &new_settings_data;
#else
    const generator_settings* new_settings = s;
#endif
    int64_t fei;
    int i;

    int64_t square_counts[GRAPHGEN_INITIATOR_SIZE2];
    make_square_counts(num_edges, &state, s, square_counts);
#ifdef GRAPHGEN_UNDIRECTED
    if (base_src == base_tgt) {
      /* Clip-and-flip for undirected graph */
      int i, j;
      for (i = 0; i < GRAPHGEN_INITIATOR_SIZE; ++i) {
        for (j = i + 1; j < GRAPHGEN_INITIATOR_SIZE; ++j) {
          square_counts[i * GRAPHGEN_INITIATOR_SIZE + j] += square_counts[j * GRAPHGEN_INITIATOR_SIZE + i];
          square_counts[j * GRAPHGEN_INITIATOR_SIZE + i] = 0;
        }
      }
    }
#endif /* GRAPHGEN_UNDIRECTED */
    new_nverts = nverts / GRAPHGEN_INITIATOR_SIZE;
    fei = first_edge_index;
    for (i = 0; i < GRAPHGEN_INITIATOR_SIZE2; ++i) {
      if (square_counts[i] != 0) {
        int lhs_pos = (fei < my_first_edge) ? -1 : (fei < my_last_edge) ? 0 : 1;
        int rhs_pos = ((fei + square_counts[i]) < my_first_edge) ? -1 : ((fei + square_counts[i]) < my_last_edge) ? 0 : 1;
        if (lhs_pos == 0 || rhs_pos == 0 || lhs_pos != rhs_pos) {
#ifdef GRAPHGEN_MODIFY_PARAMS_AT_EACH_LEVEL
          new_settings_data = *s;
          alter_params(&new_settings_data, trans + 1, &subpart_states[i]);
#endif
          generate_kronecker_internal(
            orig_state,
            fei,
            square_counts[i],
            new_nverts,
            new_settings,
            base_src + new_nverts * (i / GRAPHGEN_INITIATOR_SIZE),
            base_tgt + new_nverts * (i % GRAPHGEN_INITIATOR_SIZE));
        }
        fei += square_counts[i];
      }
    }
  }
}

int64_t compute_edge_array_size(
       int rank, int size,
       int64_t M) {
  int64_t rankc = (int64_t)(rank);
  int64_t sizec = (int64_t)(size);
  int64_t my_start_edge = rankc * (M / sizec) + (rankc < (M % sizec) ? rankc : (M % sizec));
  int64_t my_end_edge = (rankc + 1) * (M / sizec) + (rankc + 1 < (M % sizec) ? rankc + 1 : (M % sizec));
  return my_end_edge - my_start_edge;
}

void generate_kronecker(
       int rank, int size,
       const uint_fast32_t seed[5] /* All values in [0, 2^31 - 1), not all zero */,
       int logN /* In base GRAPHGEN_INITIATOR_SIZE */,
       int64_t M,
       const double initiator[GRAPHGEN_INITIATOR_SIZE2],
#ifdef GRAPHGEN_KEEP_MULTIPLICITIES
       generated_edge* const edges /* Size >= compute_edge_array_size(rank, size, M), must be zero-initialized */
#else
       int64_t* const edges /* Size >= 2 * compute_edge_array_size(rank, size, M) */
#endif
) {

  mrg_state state;
  int64_t my_start_edge = rank * (M / size) + (rank < (M % size) ? rank : (M % size));
  int64_t my_end_edge = (rank + 1) * (M / size) + (rank + 1 < (M % size) ? rank + 1 : (M % size));
  unsigned int i;
  generator_settings settings_data;

  mrg_seed(&state, seed);

  for (i = 0; i < GRAPHGEN_INITIATOR_SIZE2; ++i) {
    settings_data.initiator[i] = initiator[i];
  }
  settings_data.my_first_edge = my_start_edge;
  settings_data.my_last_edge = my_end_edge;
  settings_data.total_nverts = (int64_t)pow(GRAPHGEN_INITIATOR_SIZE, logN);
  settings_data.out = edges;

  generate_kronecker_internal(
    &state,
    0,
    M,
    (int64_t)pow(GRAPHGEN_INITIATOR_SIZE, logN),
    &settings_data,
    0,
    0);
}
