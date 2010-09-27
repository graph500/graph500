/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#if !defined(XALLOC_HEADER_)
#define XALLOC_HEADER_

void * xmalloc (size_t);
void * xmalloc_large (size_t);
void xfree_large (void *);
void * xmalloc_large_ext (size_t);

/*
void mark_large_unused (void *);
void mark_large_willuse (void *);
*/

#endif /* XALLOC_HEADER_ */
