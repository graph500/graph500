/* Copyright (C) 2010-2011 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */
/*           Anton Korzh                                                   */

#define ulong_bits 64
#define ulong_mask &63
#define ulong_shift >>6
#define SET_VISITED(v) do {visited[VERTEX_LOCAL((v)) ulong_shift] |= (1UL << (VERTEX_LOCAL((v)) ulong_mask));} while (0)
#define SET_VISITEDLOC(v) do {visited[(v) ulong_shift] |= (1ULL << ((v) ulong_mask));} while (0)
#define TEST_VISITED(v) ((visited[VERTEX_LOCAL((v)) ulong_shift] & (1UL << (VERTEX_LOCAL((v)) ulong_mask))) != 0)
#define TEST_VISITEDLOC(v) ((visited[(v) ulong_shift] & (1ULL << ((v) ulong_mask))) != 0)
#define CLEAN_VISITED()  memset(visited,0,visited_size*sizeof(unsigned long));
