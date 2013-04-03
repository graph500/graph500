#if !defined(GLOBALS_HEADER_)
#define GLOBALS_HEADER_

#if !defined(SCALE_MAX)
#define SCALE_MAX 40
#elif SCALE_MAX > 53
#error "Not using sufficient randomness in root sampling for SCALE_MAX > 53."
#endif

#if !defined(SCALE_DEFAULT)
#define SCALE_DEFAULT 12
#endif

#if !defined(EF_DEFAULT)
#define EF_DEFAULT 16
#endif

#if !defined(MAXWEIGHT_DEFAULT)
#define MAXWEIGHT_DEFAULT 255
#endif

#if !defined(NROOT_MAX)
#define NROOT_MAX 32
#endif

#if !defined(NROOT_DEFAULT)
#define NROOT_DEFAULT 16
#endif

#if !defined(A_DEFAULT)
#define A_DEFAULT 0.55
#endif

#if !defined(B_DEFAULT)
#define B_DEFAULT 0.1
#endif

#if !defined(NOISEFACT_DEFAULT)
#define NOISEFACT_DEFAULT 0.1
#endif

#if !defined(IN_GLOBALS_C)
extern const int SCALE, EF, NROOT, MAXWEIGHT;
extern const int64_t NV, NE, Z, Zinv;
extern const float A, B, NOISEFACT;
#endif

void init_globals (int, int, int, int, float, float, float);

#endif /* GLOBALS_HEADER_ */
