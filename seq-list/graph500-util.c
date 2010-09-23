/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <errno.h>

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

#include <alloca.h> /* Portable enough. */
#include <unistd.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#if !defined(MAP_HUGETLB)
#define MAP_HUGETLB 0
#endif
#if !defined(MAP_POPULATE)
#define MAP_POPULATE 0
#endif

#include "graph500-util.h"

#define NSTAT 9
#define PRINT_STATS(lbl, israte)					\
  do {									\
    printf ("min_%s: %20.17e\n", lbl, stats[0]);			\
    printf ("firstquartile_%s: %20.17e\n", lbl, stats[1]);		\
    printf ("median_%s: %20.17e\n", lbl, stats[2]);			\
    printf ("thirdquartile_%s: %20.17e\n", lbl, stats[3]);		\
    printf ("max_%s: %20.17e\n", lbl, stats[4]);			\
    if (!israte) {							\
      printf ("mean_%s: %20.17e\n", lbl, stats[5]);			\
      printf ("stddev_%s: %20.17e\n", lbl, stats[6]);			\
    } else {								\
      printf ("harmonic_mean_%s: %20.17e\n", lbl, stats[7]);		\
      printf ("harmonic_stddev_%s: %20.17e\n", lbl, stats[8]);	\
    }									\
  } while (0)

void
output_results (const int64_t SCALE, int64_t nvtx_scale, int64_t edgefactor,
		const double A, const double B, const double C, const double D,
		const double construction_time,
		const int NBFS, const double *bfs_time, const int64_t *bfs_nedge)
{
  int k;
  double *tm;
  double *stats;

  tm = alloca (NBFS * sizeof (*tm));
  stats = alloca (NSTAT * sizeof (*stats));
  if (!tm || !stats) {
    perror ("Error allocating within final statistics calculation.");
    abort ();
  }

  printf ("SCALE: %" PRId64"\nnvtx: %" PRId64 "\nedgefactor: %" PRId64 "\n",
	  SCALE, nvtx_scale, edgefactor);
  printf ("A: %20.17e\nB: %20.17e\nC: %20.17e\nD: %20.17e\n", A, B, C, D);
  printf ("construction_time: %20.17e\n", construction_time);
  printf ("nbfs: %d\n", NBFS);

  memcpy (tm, bfs_time, NBFS*sizeof(tm[0]));
  statistics (stats, tm, NBFS);
  PRINT_STATS("time", 0);

  for (k = 0; k < NBFS; ++k)
    tm[k] = bfs_nedge[k];
  statistics (stats, tm, NBFS);
  PRINT_STATS("nedge", 0);

  for (k = 0; k < NBFS; ++k)
    tm[k] = bfs_nedge[k] / bfs_time[k];
  statistics (stats, tm, NBFS);
  PRINT_STATS("TEPS", 1);
}

static int
dcmp (const void *a, const void *b)
{
  const double da = *(const double*)a;
  const double db = *(const double*)b;
  if (da > db) return 1;
  if (db > da) return -1;
  if (da == db) return 0;
  fprintf (stderr, "No NaNs permitted in output.\n");
  abort ();
}

void
statistics (double *out, double *data, int64_t n)
{
  long double s, mean;
  double t;
  int k;

  /* Quartiles */
  qsort (data, n, sizeof (*data), dcmp);
  out[0] = data[0];
  t = (n+1) / 4.0;
  k = (int) t;
  if (t == k)
    out[1] = data[k];
  else
    out[1] = 3*(data[k]/4.0) + data[k+1]/4.0;
  t = (n+1) / 2.0;
  k = (int) t;
  if (t == k)
    out[2] = data[k];
  else
    out[2] = data[k]/2.0 + data[k+1]/2.0;
  t = 3*((n+1) / 4.0);
  k = (int) t;
  if (t == k)
    out[3] = data[k];
  else
    out[3] = data[k]/4.0 + 3*(data[k+1]/4.0);
  out[4] = data[n-1];

  s = data[n-1];
  for (k = n-1; k > 0; --k)
    s += data[k-1];
  mean = s/n;
  out[5] = mean;
  s = data[n-1] - mean;
  s *= s;
  for (k = n-1; k > 0; --k) {
    long double tmp = data[k-1] - mean;
    s += tmp * tmp;
  }
  out[6] = sqrt (s/(n-1));

  s = (data[0]? 1.0L/data[0] : 0);
  for (k = 1; k < n; ++k)
    s += (data[k]? 1.0L/data[k] : 0);
  out[7] = n/s;
  mean = s/n;

  /*
    Nilan Norris, The Standard Errors of the Geometric and Harmonic
    Means and Their Application to Index Numbers, 1940.
    http://www.jstor.org/stable/2235723
  */
  s = (data[0]? 1.0L/data[0] : 0) - mean;
  s *= s;
  for (k = 1; k < n; ++k) {
    long double tmp = (data[k]? 1.0L/data[k] : 0) - mean;
    s += tmp * tmp;
  }
  s = (sqrt (s)/(n-1)) * out[7] * out[7];
  out[8] = s;
}

void *
xmalloc (size_t sz)
{
  void *out;
  if (!(out = malloc (sz))) {
    perror ("malloc failed");
    abort ();
  }
  return out;
}

void *
xmalloc_large (size_t sz)
{
#if defined(USE_MMAP_LARGE)
  void *out;
  out = mmap (NULL, sz, PROT_READ|PROT_WRITE,
	      MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_POPULATE, 0, 0);
  if (!out) {
    perror ("mmap failed");
    abort ();
  }
  return out;
#else
  return xmalloc (sz);
#endif
}

void
xfree_large (void *p, size_t sz)
{
#if defined(USE_MMAP_LARGE)
  munmap (p, sz);
#else
  free (p);
#endif
}

static int ext_fd = -1;
static void *ext_p = NULL;
static size_t ext_sz;
static void (*old_abort_handler)(int);

static void
exit_handler (void)
{
  if (ext_p)
    munmap (ext_p, ext_sz);
  if (ext_fd >= 0)
    close (ext_fd);
  ext_fd = -1;
}

static void
abort_handler (int passthrough)
{
  exit_handler ();
  old_abort_handler (passthrough);
}

void *
xmalloc_large_ext (size_t sz)
{
  char extname[] = "/tmp/graph500-IJ-XXXXXX";
  void *out;
  int fd;

  if (ext_fd >= 0) {
    fprintf (stderr, "Cannot use xmalloc_large_ext more than once.\n");
    abort ();
  }

  fd = mkstemp (extname);
  if (fd < 0) {
    perror ("xmalloc_large_ext failed to make a file");
    abort ();
  }
  if (unlink (extname)) {
    perror ("UNLINK FAILED!");
    goto errout;
  }

  if (lseek (fd, sz - sizeof(fd), SEEK_SET) < 0) {
    perror ("lseek failed");
    goto errout;
  }
  if (write (fd, &fd, sizeof(fd)) != sizeof (fd)) {
    perror ("resizing write failed");
    goto errout;
  }

  out = mmap (NULL, sz, PROT_READ|PROT_WRITE,
	      MAP_PRIVATE|MAP_HUGETLB|MAP_POPULATE, fd, 0);
  if (!out) {
    perror ("mmap failed");
    goto errout;
  }

  if (atexit (exit_handler)) {
    perror ("failed to install exit handler");
    goto errout;
  }

  old_abort_handler = signal (SIGABRT, abort_handler);
  if (SIG_ERR == old_abort_handler) {
    perror ("failed to install cleanup handler");
    goto errout;
  }

  ext_p = out;
  ext_fd = fd;
  ext_sz = sz;
  return out;

 errout:
  if (fd >= 0) close (fd);
  abort ();
}

void
xfree_large_ext (void)
{
  exit_handler ();
}

void
mark_ext_unused (void)
{
  if (ext_p)
    posix_madvise (ext_p, ext_sz, POSIX_MADV_DONTNEED);
}

void
mark_ext_willuse (void)
{
  if (ext_p) {
    posix_madvise (ext_p, ext_sz, POSIX_MADV_WILLNEED);
    posix_madvise (ext_p, ext_sz, POSIX_MADV_SEQUENTIAL);
  }
}

#if defined(CLOCK_REALTIME)
#define TICTOC_CLOCK CLOCK_REALTIME
#define TICTOC_CLOCK_NAME "CLOCK_REALTIME"
#else
#error "Failed to find a timing clock."
#endif

static clockid_t tictoc_clock = TICTOC_CLOCK;
static struct timespec tic_ts;

void
init_tictoc (void)
{
#if _XOPEN_SOURCE >= 600
  int err;
  err = clock_getcpuclockid (0, &tictoc_clock);
  if (err >= 0) return;
  fprintf (stderr, "Unable to find CPU clock, falling back to "
	   TICTOC_CLOCK_NAME "\n");
#else
  /* Just use the default... */
#endif
}

void
tic (void)
{
  clock_gettime (TICTOC_CLOCK, &tic_ts);
}

double
toc (void)
{
  struct timespec ts;
  double out;
  clock_gettime (TICTOC_CLOCK, &ts);

  out = (ts.tv_nsec - (double)tic_ts.tv_nsec) * 1.0e-9;
  out += (ts.tv_sec - (double)tic_ts.tv_sec);
  return out;
}
