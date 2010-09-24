#if !defined(COMPAT_HEADER_)
#define COMPAT_HEADER_

#if __STDC_VERSION__ >= 199901L
#include <inttypes.h>
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
#else
#define OMP(x)
#endif

#if defined(__MTA__)
#define MTA(x) _Pragma(x)
#else
#define MTA(x)
#endif

#endif /* COMPAT_HEADER_ */
