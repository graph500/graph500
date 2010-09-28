/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#if !defined(PRNG_HEADER_)
#define PRNG_HEADER_

/** Initialze the PRNG, called in a sequential context. */
void init_random (void);

/** Fill a double-precision vector with samples from U(0,1), parallel context. */
void random_vector (double *, int64_t);

#if !defined(DONT_USE_INTERNAL)
extern uint_fast32_t prng_seed[5];
extern void *prng_state;
#endif

#endif /* PRNG_HEADER_ */
