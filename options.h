#if !defined(OPTIONS_HEADER_)
#define OPTIONS_HEADER_

extern int VERBOSE;
extern int use_RMAT;

#define A_PARAM 0.57
#define B_PARAM 0.19
#define C_PARAM 0.19
/* Hence D = 0.05. */

extern double A, B, C, D;

#define NBFS_max 64
extern int NBFS;

#define default_SCALE ((int64_t)14)
#define default_edgefactor ((int64_t)16)

extern int64_t SCALE;
extern int64_t edgefactor;

void get_options (int argc, char **argv);

#endif /* OPTIONS_HEADER_ */
