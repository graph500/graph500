#if !defined(COMPAT_HEADER_)
#define COMPAT_HEADER_

#define _FILE_OFFSET_BITS 64
#define _THREAD_SAFE
#define _XOPEN_SOURCE 600
#define _XOPEN_SOURCE_EXTENDED
#define _SVID_SOURCE

#if __STDC_VERSION__ >= 199901L
#include <inttypes.h>
#elif defined(__MTA__)
#define PRId64 "ld"
#define SCNd64 "ld"
#else
#warning "Defining long as int64_t."
typedef long int64_t;
#define PRId64 "ld"
#define SCNd64 "ld"
#if !defined(restrict)
#define restrict
#endif
#endif

#if defined(_OPENMP)
#define OMP(x) _Pragma(x)
#include <omp.h>
#else
#define OMP(x)
static int omp_get_num_threads (void) { return 1; }
static int omp_get_thread_num (void) { return 0; }
#endif

#if defined(__MTA__)
#define MTA(x) _Pragma(x)
#else
#define MTA(x)
#endif

#endif /* COMPAT_HEADER_ */
