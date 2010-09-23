/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef MOD_ARITH_H
#define MOD_ARITH_H

/* Various modular arithmetic operations for modulus 2^31-1 (0x7FFFFFFF).
 * These may need to be tweaked to get acceptable performance on some platforms
 * (especially ones without conditional moves). */

/* This code is now just a dispatcher that chooses the right header file to use
 * per-platform. */

/* FIXME: fill this in automatically */
#ifndef FAST_64BIT_ARITHMETIC
#define FAST_64BIT_ARITHMETIC
#endif

#ifdef __MTA__
#include "mod_arith_xmt.h"
#else
#ifdef FAST_64BIT_ARITHMETIC
#include "mod_arith_64bit.h"
#else
#include "mod_arith_32bit.h"
#endif
#endif

#endif /* MOD_ARITH_H */
