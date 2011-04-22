/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#if !defined(PRNG_HEADER_)
#define PRNG_HEADER_

/** Initialze the PRNG, called in a sequential context. */
void init_random (void);

extern uint64_t userseed;
extern uint_fast32_t prng_seed[5];
extern void *prng_state;

#endif /* PRNG_HEADER_ */
