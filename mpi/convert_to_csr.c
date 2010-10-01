/* Copyright (C) 2010 The Trustees of Indiana University.                  */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include "common.h"
#include <mpi.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* STINGER-like group of edges, forming a linked list of pages. */
typedef struct edge_page {
  int64_t targets[16]; /* Unused elements filled with -1 */
  struct edge_page* next; /* NULL if no more */
} edge_page;

static inline edge_page* new_edge_page(void) {
  edge_page* ep = (edge_page*)xmalloc(sizeof(edge_page));
  int i;
  ep->next = NULL;
  for (i = 0; i < 16; ++i) ep->targets[i] = -1;
  return ep;
}

static inline void delete_edge_page(edge_page* ep) {
  if (!ep) return;
  delete_edge_page(ep->next);
  free(ep);
}

typedef struct adjacency_list {
  size_t nvertices; /* User-visible number of vertices */
  size_t nvertices_allocated; /* Actual allocated size of data */
  edge_page** data; /* Indexed by vertex */
} adjacency_list;

static void grow_adj_list(adjacency_list* al, size_t min_nvertices) {
  if (min_nvertices <= al->nvertices) return;
  while (min_nvertices > al->nvertices_allocated) {
    al->nvertices_allocated = (al->nvertices_allocated == 0) ? 16 : (al->nvertices_allocated * 2);
    al->data = (edge_page**)xrealloc(al->data, al->nvertices_allocated * sizeof(edge_page*));
  }
  size_t i;
  for (i = al->nvertices; i < min_nvertices; ++i) {
    al->data[i] = NULL;
  }
  al->nvertices = min_nvertices;
}

static void add_adj_list_edge(adjacency_list* al, size_t src, int64_t tgt) {
  grow_adj_list(al, src + 1);
  edge_page** p = al->data + src;
  /* Each page is filled before we allocate another one, so we only need to
   * check the last one in the chain. */
  while (*p && (*p)->next) {p = &((*p)->next);}
  if (*p) {
    assert (!(*p)->next);
    int i;
    for (i = 0; i < 16; ++i) {
      if ((*p)->targets[i] == -1) {
        (*p)->targets[i] = tgt;
        return;
      }
    }
    p = &((*p)->next);
    assert (!*p);
  }
  assert (!*p);
  *p = new_edge_page();
  (*p)->targets[0] = tgt;
}

static void clear_adj_list(adjacency_list* al) {
  size_t i;
  for (i = 0; i < al->nvertices; ++i) delete_edge_page(al->data[i]);
  free(al->data);
  al->data = NULL;
  al->nvertices = al->nvertices_allocated = 0;
}

void convert_graph_to_csr(const int64_t nedges, const int64_t* const edges, csr_graph* const g) {
  adjacency_list adj_list = {0, 0, NULL}; /* Adjacency list being built up with
                                           * received data */
  {
    /* Redistribute each input undirected edge (a, b) to create two directed
     * copies: (a -> b) on VERTEX_OWNER(a) and (b -> a) on VERTEX_OWNER(b)
     * [except for self-loops, of which only one copy is kept]. */
    const size_t edge_buffer_size = (1 << 27) / (2 * sizeof(int64_t)) / size; /* 128 MiB */
    /* Note that these buffers are edge pairs (src, tgt), where both elements
     * are global vertex indexes. */
    int64_t* recvbuf = (int64_t*)xmalloc(edge_buffer_size * 2 * sizeof(int64_t));
    MPI_Request recvreq;
    int recvreq_active = 0;
    int64_t* coalescing_buf = (int64_t*)xmalloc(size * edge_buffer_size * 2 * sizeof(int64_t));
    size_t* coalescing_counts = (size_t*)xcalloc(size, sizeof(size_t)); /* Uses zero-init */
    MPI_Request* sendreqs = (MPI_Request*)xmalloc(size * sizeof(MPI_Request));
    int* sendreqs_active = (int*)xcalloc(size, sizeof(int)); /* Uses zero-init */
    int num_sendreqs_active = 0;
    int num_senders_done = 1; /* Number of ranks that have said that they will
                               * not send to me again; I will never send to
                               * myself at all (see test in SEND). */

    /* The heavy use of macros here is to create the equivalent of nested
     * functions with proper access to variables from the enclosing scope. */

#define PROCESS_EDGE(src, tgt) \
    /* Do the handling for one received edge. */ \
    do { \
      assert (VERTEX_OWNER((src)) == rank); \
      add_adj_list_edge(&adj_list, VERTEX_LOCAL((src)), (tgt)); \
    } while (0)

#define START_IRECV \
    /* Start/restart the receive operation to wait for blocks of edges, if
     * needed. */ \
    do { \
      if (num_senders_done < size) { \
        MPI_Irecv(recvbuf, edge_buffer_size * 2, INT64_T_MPI_TYPE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recvreq); \
        recvreq_active = 1; \
      } \
    } while (0)

#define PROCESS_REQS \
    /* Handle all outstanding MPI requests and progress the MPI implementation.
     * */ \
    do { \
      int flag; \
      /* Test receive request. */ \
      while (recvreq_active) { \
        MPI_Status st; \
        MPI_Test(&recvreq, &flag, &st); \
        if (!flag) break; \
        /* A message arrived. */ \
        recvreq_active = 0; \
        int count; \
        MPI_Get_count(&st, INT64_T_MPI_TYPE, &count); \
        count /= 2; \
        if (count == 0 /* This count is used as a flag when each sender is done */ ) { \
          ++num_senders_done; \
        } else { \
          /* Process the edges in the received message. */ \
          int c; \
          for (c = 0; c < count; ++c) { \
            PROCESS_EDGE(recvbuf[c * 2], recvbuf[c * 2 + 1]); \
          } \
        } \
        START_IRECV; \
      } \
      /* Test send requests to determine when their buffers are available for
       * reuse. */ \
      int c; \
      for (c = 0; c < size; ++c) { \
        if (sendreqs_active[c]) { \
          MPI_Test(&sendreqs[c], &flag, MPI_STATUS_IGNORE); \
          if (flag) {sendreqs_active[c] = 0; --num_sendreqs_active;} \
        } \
      } \
    } while (0)

#define SEND(src, tgt) \
    do { \
      int dest = VERTEX_OWNER((src)); \
      if (dest == rank) { \
        /* Process self-sends locally. */ \
        PROCESS_EDGE((src), (tgt)); \
      } else { \
        while (sendreqs_active[dest]) PROCESS_REQS; /* Wait for send buffer to be available */ \
        /* Push onto coalescing buffer. */ \
        size_t c = coalescing_counts[dest]; \
        coalescing_buf[dest * edge_buffer_size * 2 + c * 2] = (src); \
        coalescing_buf[dest * edge_buffer_size * 2 + c * 2 + 1] = (tgt); \
        ++coalescing_counts[dest]; \
        /* Send if the buffer is full. */ \
        if (coalescing_counts[dest] == edge_buffer_size) { \
          FLUSH_COALESCING_BUFFER(dest); \
        } \
      } \
    } while (0)

#define FLUSH_COALESCING_BUFFER(dest) \
    do { \
      while (sendreqs_active[(dest)]) PROCESS_REQS; /* Wait for previous sends to finish */ \
      /* Ssend plus only having one request to a given destination active at a time should act as flow control. */ \
      MPI_Issend(coalescing_buf + (dest) * edge_buffer_size * 2, coalescing_counts[(dest)] * 2, INT64_T_MPI_TYPE, (dest), 0, MPI_COMM_WORLD, &sendreqs[(dest)]); \
      sendreqs_active[(dest)] = 1; \
      ++num_sendreqs_active; \
      /* Clear the buffer for the next user. */ \
      coalescing_counts[(dest)] = 0; \
    } while (0)

    START_IRECV;
    size_t i;
    for (i = 0; i < (size_t)nedges; ++i) {
      if ((i % (1 << 16)) == 0) PROCESS_REQS;
      if (edges[i * 2 + 0] == -1 || edges[i * 2 + 1] == -1) {
        continue;
      }
      SEND(edges[i * 2 + 0], edges[i * 2 + 1]);
      if (edges[i * 2 + 0] != edges[i * 2 + 1]) {
        /* Only send reverse for non-self-loops. */
        SEND(edges[i * 2 + 1], edges[i * 2 + 0]);
      }
    }
    int offset;
    for (offset = 1; offset < size; ++offset) {
      int dest = MOD_SIZE(rank + offset);
      if (coalescing_counts[dest] != 0) {
        /* Send actual data, if any. */
        FLUSH_COALESCING_BUFFER(dest);
      }
      /* Send empty message to indicate that we won't send anything else to
       * this rank (takes advantage of MPI non-overtaking rules). */
      FLUSH_COALESCING_BUFFER(dest);
    }
    while (num_senders_done < size || num_sendreqs_active > 0) PROCESS_REQS;
    free(recvbuf);
    free(coalescing_buf);
    free(coalescing_counts);
    free(sendreqs);
    free(sendreqs_active);

#undef PROCESS_REQS
#undef PROCESS_EDGE
#undef FLUSH_COALESCING_BUFFER
#undef SEND
#undef START_IRECV

  }

  /* Compute global number of vertices and count the degrees of the local
   * vertices. */
  int64_t nverts_known = 0; /* We only count vertices touched by at least one
                             * edge, and because of edge doubling each vertex
                             * incident to an edge must be the target of some
                             * copy of that edge. */
  size_t nlocalverts_orig = adj_list.nvertices;
  size_t* degrees = (size_t*)xcalloc(nlocalverts_orig, sizeof(size_t)); /* Uses zero-init */
  size_t i, j;
  for (i = 0; i < nlocalverts_orig; ++i) {
    size_t deg = 0;
    edge_page* p;
    for (p = adj_list.data[i]; p; p = p->next) {
      for (j = 0; j < 16; ++j) {
        if (p->targets[j] != -1) {
          ++deg;
          if (p->targets[j] >= nverts_known) nverts_known = p->targets[j] + 1;
        }
      }
    }
    degrees[i] = deg;
  }
  int64_t nglobalverts = 0;
  MPI_Allreduce(&nverts_known, &nglobalverts, 1, INT64_T_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD);
  g->nglobalverts = nglobalverts;
  /* Compute the final number of local vertices based on the global maximum
   * vertex number. */
  size_t nlocalverts = VERTEX_LOCAL(nglobalverts + size - 1 - rank);
  g->nlocalverts = nlocalverts;
  grow_adj_list(&adj_list, nlocalverts);

  /* Build CSR data structure. */
  size_t *rowstarts = (size_t*)xmalloc((nlocalverts + 1) * sizeof(size_t));
  g->rowstarts = rowstarts;
  /* Compute offset to start of each row. */
  rowstarts[0] = 0;
  for (i = 0; i < nlocalverts; ++i) {
    rowstarts[i + 1] = rowstarts[i] + (i >= nlocalverts_orig ? 0 : degrees[i]);
  }
  size_t nlocaledges = rowstarts[nlocalverts];
  g->nlocaledges = nlocaledges;
  int64_t* column = (int64_t*)xmalloc(nlocaledges * sizeof(int64_t));
  g->column = column;
  /* Append outgoing edges for each vertex to the column array, in order. */
  for (i = 0; i < nlocalverts; ++i) {
    edge_page* p;
    int offset = 0;
    for (p = adj_list.data[i]; p; p = p->next, offset += 16) {
      size_t deg = (i >= nlocalverts_orig ? 0 : degrees[i]);
      size_t nelts = (deg - offset > 16) ? 16 : deg - offset;
      memcpy(column + rowstarts[i] + offset,
             p->targets,
             nelts * sizeof(int64_t));
    }
  }
  free(degrees); degrees = NULL;
  clear_adj_list(&adj_list);
}
