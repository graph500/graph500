#if !defined(SORTS_HEADER_)
#define SORTS_HEADER_

void insertion_i64 (int64_t * restrict, const size_t);
void shellsort_i64 (int64_t * restrict, const size_t);
void introsort_i64 (int64_t * restrict, const size_t);

void insertion_both_i64 (int64_t * restrict, int16_t * restrict, const size_t);
void shellsort_both_i64 (int64_t * restrict, int16_t * restrict, const size_t);
void introsort_both_i64 (int64_t * restrict, int16_t * restrict, const size_t);

#endif /* SORTS_HEADER_ */
