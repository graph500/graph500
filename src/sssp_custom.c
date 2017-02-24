// Stub for custom SSSP implementation

#include "aml.h"
#include "common.h"
#include "csr_reference.h"
#include "bitmap_reference.h"

extern oned_csr_graph g;
extern int64_t* column,*pred_glob,visited_size;
extern unsigned long * visited;

#ifdef SSSP

//user provided function to be called several times to implement kernel 3: single source shortest path
//function has filled with -1 and 1.0 pred and dist arrays
//at exit dist should have shortest distance to root for each vertice otherwise -1
//pred array should point to next vertie in shortest path or -1 if vertex unreachable
//pred[VERTEX_LOCAL(root)] should be root and dist should be 0.0

void run_sssp(int64_t root,int64_t* pred,float *dist) {
	pred_glob=pred;
	//user code for SSSP
}

//user provided function to prefill dist array with whatever value
void clean_shortest(float* dist) {
	int i;
	for(i=0;i<g.nlocalverts;i++) dist[i]=-1.0;
}
#endif
