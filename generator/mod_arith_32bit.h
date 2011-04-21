/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef MOD_ARITH_32BIT_H
#define MOD_ARITH_32BIT_H

#include <stdint.h>
#include <assert.h>

/* Various modular arithmetic operations for modulus 2^31-1 (0x7FFFFFFF).
 * These may need to be tweaked to get acceptable performance on some platforms
 * (especially ones without conditional moves). */

static inline uint_fast32_t mod_add(uint_fast32_t a, uint_fast32_t b) {
  uint_fast32_t x;
  assert (a <= 0x7FFFFFFE);
  assert (b <= 0x7FFFFFFE);
#if 0
  return (a + b) % 0x7FFFFFFF;
#else
  x = a + b; /* x <= 0xFFFFFFFC */
  x = (x >= 0x7FFFFFFF) ? (x - 0x7FFFFFFF) : x;
  return x;
#endif
}

static inline uint_fast32_t mod_mul(uint_fast32_t a, uint_fast32_t b) {
  uint_fast64_t temp;
  uint_fast32_t temp2;
  assert (a <= 0x7FFFFFFE);
  assert (b <= 0x7FFFFFFE);
#if 0
  return (uint_fast32_t)((uint_fast64_t)a * b % 0x7FFFFFFF);
#else
  temp = (uint_fast64_t)a * b; /* temp <= 0x3FFFFFFE00000004 */
  temp2 = (uint_fast32_t)(temp & 0x7FFFFFFF) + (uint_fast32_t)(temp >> 31); /* temp2 <= 0xFFFFFFFB */
  return (temp2 >= 0x7FFFFFFF) ? (temp2 - 0x7FFFFFFF) : temp2;
#endif
}

static inline uint_fast32_t mod_mac(uint_fast32_t sum, uint_fast32_t a, uint_fast32_t b) {
  uint_fast64_t temp;
  uint_fast32_t temp2;
  assert (sum <= 0x7FFFFFFE);
  assert (a <= 0x7FFFFFFE);
  assert (b <= 0x7FFFFFFE);
#if 0
  return (uint_fast32_t)(((uint_fast64_t)a * b + sum) % 0x7FFFFFFF);
#else
  temp = (uint_fast64_t)a * b + sum; /* temp <= 0x3FFFFFFE80000002 */
  temp2 = (uint_fast32_t)(temp & 0x7FFFFFFF) + (uint_fast32_t)(temp >> 31); /* temp2 <= 0xFFFFFFFC */
  return (temp2 >= 0x7FFFFFFF) ? (temp2 - 0x7FFFFFFF) : temp2;
#endif
}

static inline uint_fast32_t mod_mac2(uint_fast32_t sum, uint_fast32_t a, uint_fast32_t b, uint_fast32_t c, uint_fast32_t d) {
  assert (sum <= 0x7FFFFFFE);
  assert (a <= 0x7FFFFFFE);
  assert (b <= 0x7FFFFFFE);
  assert (c <= 0x7FFFFFFE);
  assert (d <= 0x7FFFFFFE);
  return mod_mac(mod_mac(sum, a, b), c, d);
}

static inline uint_fast32_t mod_mac3(uint_fast32_t sum, uint_fast32_t a, uint_fast32_t b, uint_fast32_t c, uint_fast32_t d, uint_fast32_t e, uint_fast32_t f) {
  assert (sum <= 0x7FFFFFFE);
  assert (a <= 0x7FFFFFFE);
  assert (b <= 0x7FFFFFFE);
  assert (c <= 0x7FFFFFFE);
  assert (d <= 0x7FFFFFFE);
  assert (e <= 0x7FFFFFFE);
  assert (f <= 0x7FFFFFFE);
  return mod_mac2(mod_mac(sum, a, b), c, d, e, f);
}

static inline uint_fast32_t mod_mac4(uint_fast32_t sum, uint_fast32_t a, uint_fast32_t b, uint_fast32_t c, uint_fast32_t d, uint_fast32_t e, uint_fast32_t f, uint_fast32_t g, uint_fast32_t h) {
  assert (sum <= 0x7FFFFFFE);
  assert (a <= 0x7FFFFFFE);
  assert (b <= 0x7FFFFFFE);
  assert (c <= 0x7FFFFFFE);
  assert (d <= 0x7FFFFFFE);
  assert (e <= 0x7FFFFFFE);
  assert (f <= 0x7FFFFFFE);
  assert (g <= 0x7FFFFFFE);
  assert (h <= 0x7FFFFFFE);
  return mod_mac2(mod_mac2(sum, a, b, c, d), e, f, g, h);
}

/* The two constants x and y are special cases because they are easier to
 * multiply by on 32-bit systems.  They are used as multipliers in the random
 * number generator.  The techniques for fast multiplication by these
 * particular values are from:
 *
 * Pierre L'Ecuyer, Francois Blouin, and Raymond Couture. 1993. A search
 * for good multiple recursive random number generators. ACM Trans. Model.
 * Comput. Simul. 3, 2 (April 1993), 87-98. DOI=10.1145/169702.169698
 * http://doi.acm.org/10.1145/169702.169698
 *
 * Pierre L'Ecuyer. 1990. Random numbers for simulation. Commun. ACM 33, 10
 * (October 1990), 85-97. DOI=10.1145/84537.84555
 * http://doi.acm.org/10.1145/84537.84555
 */

inline uint_fast32_t mod_mul_x(uint_fast32_t a) {
  static const int32_t q = 20 /* UINT32_C(0x7FFFFFFF) / 107374182 */;
  static const int32_t r = 7  /* UINT32_C(0x7FFFFFFF) % 107374182 */;
  int_fast32_t result = (int_fast32_t)(a) / q;
  result = 107374182 * ((int_fast32_t)(a) - result * q) - result * r;
  result += (result < 0 ? 0x7FFFFFFF : 0);
  assert ((uint_fast32_t)(result) == mod_mul(a, 107374182));
  return (uint_fast32_t)result;
}

inline uint_fast32_t mod_mul_y(uint_fast32_t a) {
  static const int32_t q = 20554 /* UINT32_C(0x7FFFFFFF) / 104480 */;
  static const int32_t r = 1727  /* UINT32_C(0x7FFFFFFF) % 104480 */;
  int_fast32_t result = (int_fast32_t)(a) / q;
  result = 104480 * ((int_fast32_t)(a) - result * q) - result * r;
  result += (result < 0 ? 0x7FFFFFFF : 0);
  assert ((uint_fast32_t)(result) == mod_mul(a, 104480));
  return (uint_fast32_t)result;
}

inline uint_fast32_t mod_mac_y(uint_fast32_t sum, uint_fast32_t a) {
  uint_fast32_t result = mod_add(sum, mod_mul_y(a));
  assert (result == mod_mac(sum, a, 104480));
  return result;
}

#endif /* MOD_ARITH_32BIT_H */
