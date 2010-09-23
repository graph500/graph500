/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#if !defined(GRAPH500_PRNG_HEADER_)
#define GRAPH500_PRNG_HEADER_

/** Initialze the PRNG, called in a sequential context. */
void init_random (void);

/** Fill a double-precision vector with samples from U(0,1), parallel context. */
void random_vector (double *, int64_t);

#endif /* GRAPH500_PRNG_HEADER_ */
