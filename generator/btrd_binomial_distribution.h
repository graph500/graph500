/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef BTRD_BINOMIAL_H
#define BTRD_BINOMIAL_H

/* BTRD algorithm from pages 6--7 of "The Generation of Binomial Random */
/* Variates" (Wolfgang Hoermann) --                                     */
/* http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.47.8407      */

#include <math.h>
#include <stddef.h>
#include "splittable_mrg.h"

#ifdef __cplusplus
extern "C" {
#endif
size_t btrd_binomial_distribution(size_t n, double p, mrg_state* state);
#ifdef __cplusplus
}
#endif

#endif /* BTRD_BINOMIAL_H */
