#if !defined(PRNG_HEADER_)
#define PRNG_HEADER_

void init_prng (void);
int64_t scramble (int64_t);
uint8_t random_weight (int64_t);
void random_edgevals (float *, int64_t);
void sample_roots (int64_t *, int64_t, int64_t);
int32_t prng_check (void);

#endif /* PRNG_HEADER_ */
