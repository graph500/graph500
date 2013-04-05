#define _FILE_OFFSET_BITS 64
#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>

#include <errno.h>
#include <assert.h>

#include "sorts.h"

#if defined(__MTA__)
#define MTA(x) _Pragma(x)
#else
#define MTA(x)
#endif

/*
  All integer compares are in increasing order.
*/

MTA("mta expect parallel context")
void
insertion_i64 (int64_t * restrict d, const size_t N)
{
  for (size_t i = 1; i < N; ++i) {
    const int64_t key = d[i];
    size_t k;
    for (k = i; k > 0 && d[k-1] > key; --k)
      d[k] = d[k-1];
    if (i != k)
      d[k] = key;
  }
}

MTA("mta expect parallel context")
void
shellsort_i64 (int64_t * restrict d, const size_t N)
{
  /* Using Knuth's / Incerpi-Sedgewick's strides: */
  /* const size_t stride[] = { 1391376, 463792, 198768, 86961, 33936, */
  /* 				     13776, 4592, 1968, 861, 336, */
  /* 				     112, 48, 21, 7, 3, 1 }; */
  /* Using Marcin Ciura's strides: */
  const size_t stride[] = {701, 301, 132, 57, 23, 10, 4, 1};

  for (size_t k = 0; k < sizeof(stride)/sizeof(*stride); ++k) {
    const size_t h = stride[k];
    for (size_t i = h; i < N; ++i) {
      size_t j;
      int64_t t0 = d[i];
      for (j = i; j >= h && d[j-h] > t0; j -= h)
	d[j] = d[j-h];
      if (j != i)
	d[j] = t0;
    }
  }
}

#if !defined(SORT_THRESH)
#define SORT_THRESH 32
#endif

MTA("mta expect parallel context")
static inline void
thresh_shellsort_i64 (int64_t * restrict d, const size_t N)
{
  /* Using Knuth's / Incerpi-Sedgewick's strides: */
  /* const size_t stride[] = { /\* 1391376, 463792, 198768, 86961, 33936, */
  /* 			       13776, 4592, 1968, 861, 336, */
  /* 			       112, 48, *\/ 21, 7, 3, 1 }; */
  /* Using Marcin Ciura's strides: */
  const size_t stride[] = { /* 701, 301, 132, 57, */ 23, 10, 4, 1};

  for (size_t k = 0; k < sizeof(stride)/sizeof(*stride); ++k) {
    const size_t h = stride[k];
    for (size_t i = h; i < N; ++i) {
      size_t j;
      int64_t t0 = d[i];
      for (j = i; j >= h && d[j-h] > t0; j -= h)
	d[j] = d[j-h];
      if (j != i)
	d[j] = t0;
    }
  }
}

MTA("mta expect parallel context")
static int64_t
median_of_three (int64_t a, int64_t b, int64_t c)
{
  int64_t tmp;
  /* Insertion sort. */
  if (a > b) {
    tmp = a;
    a = b;
    b = tmp;
  }
  if (b > c) {
    tmp = b;
    b = c;
    c = tmp;
  }
  if (a > b) {
    tmp = a;
    a = b;
    b = tmp;
  }
  return b;
}

MTA("mta expect parallel context")
static int
floor_log2 (size_t x)
{
  if (!x) return 0;
#if defined(__GCC__)
#if SIZE_MAX < INT64_MAX
  return 31 - __builtin_clzl(x);
#else
  return 63 - __builtin_clzl(x);
#endif
#else
  int out = 0;
#if SIZE_MAX < INT64_MAX
  const uint32_t one = 1;
#else
  const uint64_t one = 1;
  if (x >= (one << 32)) { x >>= 32; out += 32; }
#endif
  if (x >= (one << 16)) { x >>= 16; out += 16; }
  if (x >= (one << 8)) { x >>=  8; out +=  8; }
  if (x >= (one << 4)) { x >>=  4; out +=  4; }
  if (x >= (one << 2)) { x >>=  2; out +=  2; }
  if (x >= (one << 1)) { out +=  1; }
  return out;
#endif
}

MTA("mta expect parallel context")
static int64_t
partition_i64 (int64_t * restrict d, int64_t i, int64_t j,
	       const int64_t piv)
{
  while (i <= j) {
    while (d[i] < piv)
      i++;
    while (d[j] > piv)
      j--;
    if (i <= j) {
      int64_t tmp0 = d[i];
      d[i] = d[j];
      d[j] = tmp0;
      i++;
      j--;
    }
  };
  return j;
}

MTA("mta expect parallel context")
static void
introsort_loop_i64 (int64_t * restrict d, const size_t start, size_t end,
		    int depth_limit)
{
  while (end - start > SORT_THRESH) {
    if (!depth_limit) {
      thresh_shellsort_i64 (&d[start], end-start); /* Yeah, should be heapsort... */
      return;
    }
    --depth_limit;
    const size_t p = partition_i64 (d, start, end-1,
				    median_of_three (d[start],
						     d[start+(end-start)/2],
						     d[end-1]));
    introsort_loop_i64 (d, p, end, depth_limit);
    end = p;
  }
}

MTA("mta expect parallel context")
void
introsort_i64 (int64_t * restrict d, const size_t N)
{
  if (N < 2) return;
  introsort_loop_i64 (d, 0, N, 2*floor_log2 (N));
  insertion_i64 (d, N);
}
