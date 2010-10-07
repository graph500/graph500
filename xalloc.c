/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#include "compat.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>
#include <time.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>

#if !defined(MAP_HUGETLB)
#define MAP_HUGETLB 0
#endif
#if !defined(MAP_POPULATE)
#define MAP_POPULATE 0
#endif
#if !defined(MAP_NOSYNC)
#define MAP_NOSYNC 0
#endif

#if 0
/* Included in the generator. */
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
#else
extern void *xmalloc (size_t);
#endif

#if defined(__MTA__)||defined(USE_MMAP_LARGE)||defined(USE_MMAP_LARGE_EXT)
#define MAX_LARGE 32
static int n_large_alloc = 0;
static struct {
  void * p;
  size_t sz;
  int fd;
} large_alloc[MAX_LARGE];

static int installed_handler = 0;
static void (*old_abort_handler)(int);

static void
exit_handler (void)
{
  int k;
  for (k = 0; k < n_large_alloc; ++k) {
    if (large_alloc[k].p)
      munmap (large_alloc[k].p, large_alloc[k].sz);
    if (large_alloc[k].fd >= 0)
      close (large_alloc[k].fd);
    large_alloc[k].p = NULL;
    large_alloc[k].fd = -1;
  }
}

static void
abort_handler (int passthrough)
{
  exit_handler ();
  if (old_abort_handler) old_abort_handler (passthrough);
}
#endif

#if !defined(MAP_ANONYMOUS) && defined(MAP_ANON)
#define MAP_ANONYMOUS MAP_ANON
#endif

void *
xmalloc_large (size_t sz)
{
#if defined(__MTA__)||defined(USE_MMAP_LARGE)
  void *out;
  int which = n_large_alloc++;
  if (n_large_alloc > MAX_LARGE) {
    fprintf (stderr, "Too many large allocations. %d %d\n", n_large_alloc, MAX_LARGE);
    --n_large_alloc;
    abort ();
  }
  large_alloc[which].p = NULL;
  large_alloc[which].fd = -1;
  out = mmap (NULL, sz, PROT_READ|PROT_WRITE,
	      MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_POPULATE, 0, 0);
  if (out == MAP_FAILED || !out) {
    perror ("mmap failed");
    abort ();
  }
  large_alloc[which].p = out;
  large_alloc[which].sz = sz;
  return out;
#else
  return xmalloc (sz);
#endif
}

void
xfree_large (void *p)
{
#if defined(__MTA__)||defined(USE_MMAP_LARGE)||defined(USE_MMAP_LARGE_EXT)
  int k, found = 0;
  for (k = 0; k < n_large_alloc; ++k) {
    if (p == large_alloc[k].p) {
      munmap (p, large_alloc[k].sz);
      large_alloc[k].p = NULL;
      if (large_alloc[k].fd >= 0) {
	close (large_alloc[k].fd);
	large_alloc[k].fd = -1;
      }
      found = 1;
      break;
    }
  }
  if (found) {
    --n_large_alloc;
    for (; k < n_large_alloc; ++k)
      large_alloc[k] = large_alloc[k+1];
  } else
    free (p);
#else
  free (p);
#endif
}

void *
xmalloc_large_ext (size_t sz)
{
#if !defined(__MTA__)&&defined(USE_MMAP_LARGE_EXT)
  char extname[PATH_MAX+1];
  char *tmppath;
  void *out;
  int fd, which;

  if (getenv ("TMPDIR"))
    tmppath = getenv ("TMPDIR");
  else if (getenv ("TEMPDIR"))
    tmppath = getenv ("TEMPDIR");
  else
    tmppath = "/tmp";

  sprintf (extname, "%s/graph500-ext-XXXXXX", tmppath);

  which = n_large_alloc++;
  if (n_large_alloc > MAX_LARGE) {
    fprintf (stderr, "Out of large allocations.\n");
    abort ();
  }
  large_alloc[which].p = 0;
  large_alloc[which].fd = -1;

  fd = mkstemp (extname);
  if (fd < 0) {
    perror ("xmalloc_large_ext failed to make a file");
    abort ();
  }
  if (unlink (extname)) {
    perror ("UNLINK FAILED!");
    goto errout;
  }

#if _XOPEN_SOURCE >= 500
  if (pwrite (fd, &fd, sizeof (fd), sz - sizeof(fd)) != sizeof (fd)) {
    perror ("resizing pwrite failed");
    goto errout;
  }
#else
  if (lseek (fd, sz - sizeof(fd), SEEK_SET) < 0) {
    perror ("lseek failed");
    goto errout;
  }
  if (write (fd, &fd, sizeof(fd)) != sizeof (fd)) {
    perror ("resizing write failed");
    goto errout;
  }
#endif
  fcntl (fd, F_SETFD, O_ASYNC);

  out = mmap (NULL, sz, PROT_READ|PROT_WRITE,
	      MAP_SHARED|MAP_POPULATE|MAP_NOSYNC, fd, 0);
  if (MAP_FAILED == out || !out) {
    perror ("mmap ext failed");
    goto errout;
  }

  if (!installed_handler) {
    installed_handler = 1;
    if (atexit (exit_handler)) {
      perror ("failed to install exit handler");
      goto errout;
    }

    old_abort_handler = signal (SIGABRT, abort_handler);
    if (SIG_ERR == old_abort_handler) {
      perror ("failed to install cleanup handler");
      goto errout;
    }
  }

  large_alloc[which].p = out;
  large_alloc[which].sz = sz;
  large_alloc[which].fd = fd;

  return out;

 errout:
  if (fd >= 0) close (fd);
  abort ();
#else
  return xmalloc_large (sz);
#endif
}

/*
void
mark_large_unused (void *p)
{
#if !defined(__MTA__)
  int k;
  for (k = 0; k < n_large_alloc; ++k)
    if (p == large_alloc[k].p)
      posix_madvise (large_alloc[k].p, large_alloc[k].sz, POSIX_MADV_DONTNEED);
#endif
}

void
mark_large_willuse (void *p)
{
#if !defined(__MTA__)
  int k;
  for (k = 0; k < n_large_alloc; ++k)
    if (p == large_alloc[k].p)
      posix_madvise (large_alloc[k].p, large_alloc[k].sz, POSIX_MADV_WILLNEED);
#endif
}
*/
