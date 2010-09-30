/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef MAKE_GRAPH_H
#define MAKE_GRAPH_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Simplified interface for users; implemented in different ways on different
 * platforms. */
void make_graph(
  /* in */ int log_numverts          /* log_2 of vertex count */,
  /* in */ int64_t desired_nedges    /* Target number of edges (actual number
                                      * will be slightly smaller) */,
  /* in */ uint64_t userseed1        /* Arbitrary 64-bit seed value */,
  /* in */ uint64_t userseed2        /* Arbitrary 64-bit seed value */,
  /* in */ const double initiator[4] /* Kronecker initiator (i.e., R-MAT a, b,
                                      * c, d) */,
  /* out */ int64_t* nedges          /* Number of generated edges */,
  /* out */ int64_t** result         /* Array of edges (each pair of elements
                                      * is a single edge); pairs with first
                                      * element -1 should be ignored;
                                      * allocated by make_graph() but must be
                                      * freed using free() by user */
);

/* PRNG interface for implementations; takes seed in same format as given by
 * users, and creates a vector of doubles in a reproducible (and
 * random-access) way. */
void make_random_numbers(
       /* in */ int64_t nvalues    /* Number of values to generate */,
       /* in */ uint64_t userseed1 /* Arbitrary 64-bit seed value */,
       /* in */ uint64_t userseed2 /* Arbitrary 64-bit seed value */,
       /* in */ int64_t position   /* Start index in random number stream */,
       /* out */ double* result    /* Returned array of values */
);

#ifdef __cplusplus
}
#endif

#endif /* MAKE_GRAPH_H */
