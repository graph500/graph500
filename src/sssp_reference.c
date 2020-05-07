/* Copyright (c) 2011-2017 Graph500 Steering Committee
   All rights reserved.
   Developed by:                Anton Korzh anton@korzh.us
                                Graph500 Steering Committee
                                http://www.graph500.org
   New code under University of Illinois/NCSA Open Source License
   see license.txt or https://opensource.org/licenses/NCSA
*/

// Graph500: Kernel 3 SSSP
// Simple parallel delta-stepping with relaxations as Active Messages

#include <string.h>
#include "aml.h"
#include "bitmap_reference.h"
#include "common.h"
#include "csr_reference.h"

#ifdef DEBUGSTATS
extern int64_t nbytes_sent, nbytes_rcvd;
#endif
// variables shared from bfs_reference
extern oned_csr_graph g;
extern int *q1, q1c;  // current bucket
extern int *q2, q2c;  // next bucket
extern int64_t *column, *glob_pred;
extern unsigned int* rowstarts;
extern unsigned long* visited;
extern int64_t visited_size;
#ifdef SSSP
// global variables as those accesed by active message handler
float* glob_dist;
float glob_maxdelta, glob_mindelta;  // range for current bucket
float* weights;
volatile int lightphase;

// Relaxation data type
typedef struct __attribute__((__packed__)) relaxmsg {
  float w;        // weight of an edge
  int dest_vloc;  // local index of destination vertex
  int src_vloc;   // local index of source vertex
} relaxmsg;

// Active message handler for relaxation
void relaxhndl(int from, void* dat, int sz) {
  relaxmsg* m = (relaxmsg*)dat;
  int vloc = m->dest_vloc;
  float w = m->w;
  float* dest_dist = &glob_dist[vloc];
  // check if relaxation is needed: either new path is shorter or vertex not
  // reached earlier
  if (*dest_dist < 0 || *dest_dist > w) {
    *dest_dist = w;                                         // update distance
    glob_pred[vloc] = VERTEX_TO_GLOBAL(from, m->src_vloc);  // update path

    if (lightphase &&
        !TEST_VISITEDLOC(vloc)) {  // Bitmap used to track if was
                                   // already relaxed with light edge
      if (w < glob_maxdelta) {     // if falls into current bucket needs further
                                   // reprocessing
        q2[q2c++] = vloc;
        SET_VISITEDLOC(vloc);
      }
    }
  }
}

// Sending relaxation active message
void send_relax(int64_t toglob, float weight, int fromloc) {
  relaxmsg m = {weight, VERTEX_LOCAL(toglob), fromloc};
  aml_send(&m, 1, sizeof(relaxmsg), VERTEX_OWNER(toglob));
}

void run_sssp(int64_t root, int64_t* pred, float* dist) {
  float delta = 0.1;
  glob_mindelta = 0.0;
  glob_maxdelta = delta;
  glob_dist = dist;
  glob_pred = pred;
  weights = g.weights;
  q1c = 0;
  q2c = 0;

  aml_register_handler(relaxhndl, 1);

  if (VERTEX_OWNER(root) == my_pe()) {
    q1[0] = VERTEX_LOCAL(root);
    q1c = 1;
    dist[VERTEX_LOCAL(root)] = 0.0;
    pred[VERTEX_LOCAL(root)] = root;
  }

  aml_barrier();

  long sum = 1;
  int64_t lastvisited = 1;
  while (sum != 0) {
#ifdef DEBUGSTATS
    double t0 = aml_time();
    nbytes_sent = 0;
#endif
    // 1. iterate over light edges: multiple phases
    while (sum != 0) {
      CLEAN_VISITED();
      lightphase = 1;
      aml_barrier();
      for (int i = 0; i < q1c; i++) {
        for (unsigned int j = rowstarts[q1[i]]; j < rowstarts[q1[i] + 1]; j++)
          if (weights[j] < delta)
            send_relax(COLUMN(j), dist[q1[i]] + weights[j], q1[i]);
      }
      aml_barrier();

      q1c = q2c;
      q2c = 0;
      int* tmp = q1;
      q1 = q2;
      q2 = tmp;
      sum = q1c;
      aml_long_allsum(&sum);
    }
    lightphase = 0;
    aml_barrier();

    // 2. iterate over S and heavy edges: one phase
    for (size_t i = 0; i < g.nlocalverts; i++)
      if (dist[i] >= glob_mindelta && dist[i] < glob_maxdelta) {
        for (unsigned int j = rowstarts[i]; j < rowstarts[i + 1]; j++)
          if (weights[j] >= delta)
            send_relax(COLUMN(j), dist[i] + weights[j], i);
      }
    aml_barrier();

    // 3. Bucket processing and checking termination condition
    glob_mindelta = glob_maxdelta;
    glob_maxdelta += delta;
    q1c = 0;
    sum = 0;
    int64_t lvlvisited = 0;
    for (size_t i = 0; i < g.nlocalverts; i++) {
      if (dist[i] >= glob_mindelta) {
        sum++;  // how many are still to be processed
        if (dist[i] < glob_maxdelta) {
          q1[q1c++] = i;  // this is lowest bucket
        }
      } else if (dist[i] != -1.0) {
        lvlvisited++;
      }
    }
    aml_long_allsum(&sum);

#ifdef DEBUGSTATS
    t0 -= aml_time();
    aml_long_allsum(&lvlvisited);
    aml_long_allsum(&nbytes_sent);
    if (!my_pe())
      printf(
          "--lvl[%1.2f..%1.2f] visited %lld (total %llu) in %5.2fs, network "
          "aggr %5.2fGb/s\n",
          glob_mindelta, glob_maxdelta, lvlvisited - lastvisited, lvlvisited,
          -t0, -(double)nbytes_sent * 8.0 / (1.e9 * t0));
    lastvisited = lvlvisited;
#endif
  }
}

void clean_shortest(float* dist) {
  for (size_t i = 0; i < g.nlocalverts; i++) dist[i] = -1.0;
}
#endif
