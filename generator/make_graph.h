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

/* Simplified interface for users; implemented in different ways on different
 * platforms. */
void make_graph(/* in */ int log_numverts   /* log_2 of vertex count */,
                /* in */ uint64_t userseed1 /* Arbitrary 64-bit seed value */,
                /* in */ uint64_t userseed2 /* Arbitrary 64-bit seed value */,
                /* out */ uint64_t* nedges  /* Number of generated edges */,
                /* out */ uint64_t** result /* Array of edges (each pair of
                                             * elements is a single edge);
                                             * pairs with first element
                                             * (uint64_t)(-1) should be
                                             * ignored; allocated by
                                             * make_graph() but must be freed
                                             * using free() by user */
);

#endif /* MAKE_GRAPH_H */
