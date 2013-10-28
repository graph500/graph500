/* Copyright 2013,  Regents of the University of California, USA. */
/* See COPYING for license. */
#ifndef BITMAP_H
#define BITMAP_H

#include "../compat.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

/*  Implemented in .h (without a .c) to allow optimizing compiler to inline
    For best performance, compile with -O3 or equivalent                     */

extern int int64_cas(int64_t* p, int64_t oldval, int64_t newval);

#define WORD_OFFSET(n) (n/64)
#define BIT_OFFSET(n) (n & 0x3f)

typedef struct {
  uint64_t *start;
  uint64_t *end;
} bitmap_t;

static inline void
bm_reset(bitmap_t* bm)
{
  OMP("omp for")
    for(uint64_t* it=bm->start; it<bm->end; it++)
      *it = 0;
}

static inline void
bm_init(bitmap_t* bm, int size)
{
  int num_longs = (size + 63) / 64;
  bm->start = (uint64_t*) malloc(sizeof(uint64_t) * num_longs);
  bm->end = bm->start + num_longs;
  bm_reset(bm);
}

static inline uint64_t
bm_get_bit(bitmap_t* bm, long pos)
{
  return bm->start[WORD_OFFSET(pos)] & (1l<<BIT_OFFSET(pos));
}

static inline long
bm_get_next_bit(bitmap_t* bm, long pos)
{
  long next = pos;
  int bit_offset = BIT_OFFSET(pos);
  uint64_t *it = bm->start + WORD_OFFSET(pos);
  uint64_t temp = (*it);
  if (bit_offset != 63) {
    temp = temp >> (bit_offset+1);
  } else {
    temp = 0;
  }
  if (!temp) {
    next = (next & 0xffffffc0);
    while (!temp) {
      it++;
      if (it >= bm->end)
        return -1;
      temp = *it;
      next += 64;
    }
  } else {
    next++;
  }
  while(!(temp&1)) {
    temp = temp >> 1;
    next++;
  }
  return next;
}

static inline void
bm_set_bit(bitmap_t* bm, long pos)
{
  bm->start[WORD_OFFSET(pos)] |= ((uint64_t) 1l<<BIT_OFFSET(pos));
}

static inline void
bm_set_bit_atomic(bitmap_t* bm, long pos)
{
  uint64_t old_val, new_val;
  uint64_t *loc = bm->start + WORD_OFFSET(pos);
  do {
    old_val = *loc;
    new_val = old_val | ((uint64_t) 1l<<BIT_OFFSET(pos));
  } while(!int64_cas((int64_t*) loc, old_val, new_val));
}

static inline void bm_swap(bitmap_t* a, bitmap_t* b)
{
  uint64_t* temp;
  temp = a->start;
  a->start = b->start;
  b->start = temp;
  temp = a->end;
  a->end = b->end;
  b->end = temp;
}

static inline void
bm_free(bitmap_t* bm)
{
  free(bm->start);
}

#endif // BITMAP_H
