/* Copyright (c) 2011-2017 Graph500 Steering Committee
   All rights reserved.
   Developed by:                Anton Korzh anton@korzh.us
                                Graph500 Steering Committee
                                http://www.graph500.org
   New code under University of Illinois/NCSA Open Source License
   see license.txt or https://opensource.org/licenses/NCSA
*/

// Graph500: Validation kernel for both BFS and SSSP
//  1. building CRS graph representation from tuple graph,
//  2. creating distances map (for BFS)
//  3. checking all edges to follow the triangle rule
//  4. checking all claimed vertices actually reached with correct distance
//  5. validating number of edges visited

#include "common.h"
#include "aml.h"
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <search.h>
#include <float.h>

#ifdef REUSE_CSR_FOR_VALIDATION
#include "csr_reference.h"
extern int64_t* column;
extern unsigned int* rowstarts;
#ifdef SSSP
extern float* weights;
#endif
#else
#define SETCOLUMN(a,b) vcolumn[a]=b
#define COLUMN(i) vcolumn[i]
#endif
int firstvalidationrun=1;
//int failedttovalidate=0;
int validatingbfs=0;

unsigned int *vdegrees,*vrowstarts;
int64_t *vcolumn;
#ifdef SSSP
float* vweights;
#endif
int64_t *globpred,nedges_traversed;
float *globdist,prevlevel;
int64_t val_errors=0;
int failedtovalidate=0;
int64_t newvisits;
int * confirmed=NULL;
int64_t maxvertex;

void frompredhndl(int from,void* data,int sz) {
	int vfrom = *(int*)data;
	int64_t predfrom = VERTEX_TO_GLOBAL(from,vfrom);
	int vloc = *(int*)(data+4);

	if(globpred[vloc] == predfrom && globdist[vloc]==FLT_MAX) globdist[vloc]=prevlevel+1.0,newvisits++;
}

void send_frompred (int vfrom,int64_t src) {
	int pe=VERTEX_OWNER(src);
	int vloc[2]={vfrom, VERTEX_LOCAL(src)};
	aml_send(&vloc,1,8,pe);
}

void vhalfedgehndl(int from,void* data,int sz)
{  vdegrees[*(int*)data]++; }

void send_half (int64_t src) {
	int pe=VERTEX_OWNER(src);
	int vloc=VERTEX_LOCAL(src);
	aml_send(&vloc,1,4,pe);
}

void vfulledgehndl(int frompe,void* data,int sz) {
	int vloc = *(int*)data;
	int64_t gtgt = *((int64_t*)(data+4));
	int next = vdegrees[vloc]++;
	SETCOLUMN(next,gtgt);
#ifdef SSSP
	vweights[next] = ((float*)data)[3];
#endif
}

void vsend_full_edge (int64_t src,int64_t tgt,float w) {
	int pe=VERTEX_OWNER(src);
	int vloc[4];
	vloc[0]=VERTEX_LOCAL(src);
	memcpy(vloc+1,&tgt,8);
#ifdef SSSP
	memcpy(vloc+3,&w,4);
	aml_send(vloc,1,16,pe);
#else
	aml_send(vloc,1,12,pe);
#endif
}

typedef struct edgedist {
	unsigned int vfrom;
	unsigned int vloc;
	int64_t predfrom;
	float distfrom;
#ifdef SSSP
	float w;
#endif
} edgedist;

void sendedgepreddist(unsigned vloc,unsigned int vedge) {
	edgedist m = {vloc,VERTEX_LOCAL(COLUMN(vedge)),globpred[vloc],globdist[vloc]
#ifdef SSSP
		,vweights[vedge]
#endif
	};
	aml_send(&m,1,sizeof(edgedist),VERTEX_OWNER(COLUMN(vedge)));
}

#define DUMPERROR(text) { printf("Validation Error: %s, edge %ld %ld weight %f pred0 %ld pred1 %ld dist0 %f dist1 %f\n",text,v0,v1,w,predv0,predv1,distv0,distv1); val_errors++; return; }

//main validation handler: tracks all edges and at delivery has both vertex preds and distances to be checked
void edgepreddisthndl(int frompe,void* data,int sz) {
	edgedist *m = (edgedist*) data;

	unsigned int v1loc = m->vloc;

	int64_t v1 = VERTEX_TO_GLOBAL(my_pe(),v1loc);
	int64_t v0 = VERTEX_TO_GLOBAL(frompe,m->vfrom);
	int64_t predv1 = globpred[v1loc];
	int64_t predv0 = m->predfrom;
#ifdef SSSP
	float w = (validatingbfs?1.0:m->w);
#else
	float w = 1.0;
#endif
	float distv0 = m->distfrom;
	float distv1 = globdist[v1loc];
	if(predv0==-1 && predv1==-1) return; else if(v0<v1) nedges_traversed++;

	if((predv0==-1 && predv1!=-1) || (predv1==-1 && predv0!=-1)) DUMPERROR("edge connecting visited and unvisited vertices");
	if(predv1==v0 && distv1 == distv0+w) confirmed[v1loc]=1; //confirm pred/dist as existing edge

	if(distv0+w < distv1 || distv1+w<distv0) DUMPERROR("triangle rule violated");
}

void makedepthmapforbfs(const size_t nlocalverts,const int64_t root,int64_t * const pred,float* dist) {

	int i,j;
	for(i=0;i<nlocalverts;i++) {
		dist[i]=FLT_MAX; //at the end there should be no FLT_MAX left
		if(pred[i]==-1) dist[i]=-1.0;
		if(pred[i]==root) dist[i]=1.0;
	}

	//fix root distance
	if(my_pe()==VERTEX_OWNER(root)) dist[VERTEX_LOCAL(root)]=0.0;

	newvisits=1;
	prevlevel=0.0;
	aml_register_handler(frompredhndl,1);

	while(newvisits!=0) {
		newvisits=0;
		prevlevel+=1.0;

		for(i=0;i<nlocalverts;i++)
			if(dist[i]==prevlevel)
				for(j=vrowstarts[i];j<vrowstarts[i+1];j++)
					send_frompred(i,COLUMN(j));
		aml_barrier();

		aml_long_allsum(&newvisits);
	}
}

int validate_result(int isbfs,const tuple_graph* const tg, const size_t nlocalverts, const int64_t root, int64_t* const pred, float* dist,int64_t *nedges_in) {
	validatingbfs=isbfs;

	//if(failedtovalidate) return 0; //failed to allocate lots of memory for validation: skipping all validation now

	size_t i,j;
	if(firstvalidationrun) {
		firstvalidationrun=0;
		confirmed = xmalloc(nlocalverts*sizeof(int));
#ifdef REUSE_CSR_FOR_VALIDATION
vrowstarts=rowstarts;
#ifdef SSSP
vweights=weights;
#endif
#else
		vdegrees=xcalloc(nlocalverts,sizeof(int));

		aml_register_handler(vhalfedgehndl,1);

		int numiters=ITERATE_TUPLE_GRAPH_BLOCK_COUNT(tg);
		// First pass : calculate degrees of each vertex
		ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize,wbuf) {
			ptrdiff_t j;
			for (j = 0; j < bufsize; ++j) {
				int64_t v0 = get_v0_from_edge(&buf[j]);
				int64_t v1 = get_v1_from_edge(&buf[j]);
				if(v0==v1) continue;
				send_half(v0);
				send_half(v1);
			}
			aml_barrier();
		} ITERATE_TUPLE_GRAPH_END;

		vrowstarts = xmalloc((nlocalverts + 1) * sizeof(int));

		vrowstarts[0] = 0;
		for (i = 0; i < nlocalverts; ++i) {
			vrowstarts[i + 1] = vrowstarts[i] + (i >= nlocalverts ? 0 : vdegrees[i]);
			vdegrees[i] = vrowstarts[i];
		}

		vcolumn = xmalloc(8*vrowstarts[nlocalverts]);
#ifdef SSSP
		vweights = xmalloc(4*vrowstarts[nlocalverts]);
#endif
		aml_register_handler(vfulledgehndl,1);
		//Second pass , actual data transfer: placing edges to its places in vcolumn
		ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize,wbuf) {
			ptrdiff_t j;
			for (j = 0; j < bufsize; ++j) {
				int64_t v0 = get_v0_from_edge(&buf[j]);
				int64_t v1 = get_v1_from_edge(&buf[j]);
				if(v0==v1) continue;
#ifdef SSSP
				vsend_full_edge(v0, v1,wbuf[j]);
				vsend_full_edge(v1, v0,wbuf[j]);
#else
				vsend_full_edge(v0, v1,1.0);
				vsend_full_edge(v1, v0,1.0);
#endif
			}
			aml_barrier();
		} ITERATE_TUPLE_GRAPH_END;

		aml_barrier();
		free(vdegrees); vdegrees=NULL;
#endif
	for (i = 0; i < nlocalverts; ++i)
		if(VERTEX_TO_GLOBAL(my_pe(),i)>maxvertex) maxvertex = VERTEX_TO_GLOBAL(my_pe(),i);
	aml_long_allmax(&maxvertex);

	} //only run at first validation

	//Actual validation here:
	globpred=pred;
	globdist=dist;
	for (i = 0; i < nlocalverts; ++i)
		confirmed[i]=0;

	for (i = 0; i < nlocalverts; ++i)
		if((pred[i]!=-1 && pred[i]<0) || pred[i]>maxvertex)
		printf("Validation Error: predecessor %ld of vertex %ld is out of range\n",pred[i],VERTEX_TO_GLOBAL(my_pe(),i)),val_errors++;
	aml_long_allsum(&val_errors);
	if(val_errors>0) return 0;

	if(validatingbfs)
	makedepthmapforbfs(nlocalverts,root,pred,dist);

	for (i = 0; i < nlocalverts; ++i) {
		if(dist[i]!=-1.0 && dist[i]<0.0)
		printf("Validation Error: distance/depth %3.2f of vertex %ld is out of range\n",dist[i],VERTEX_TO_GLOBAL(my_pe(),i)),val_errors++;
		if((pred[i]==-1 && dist[i]!=-1.0) || (pred[i]!=-1 && dist[i]==-1.0))
		printf("Validation Error: vertex %ld has inconsistent predecessor %ld and distance %f\n",VERTEX_TO_GLOBAL(my_pe(),i),pred[i],dist[i]),val_errors++;
	}

	if(my_pe()==VERTEX_OWNER(root)) {
		int vloc = VERTEX_LOCAL(root);
		if(pred[vloc]!=root || dist[vloc]!=0.0)
		printf("Validation Error: root vertex %ld has predecessor %ld and distance %f\n",root,pred[vloc],dist[vloc]),val_errors++;
		else confirmed[vloc]=1;
	}

	aml_register_handler(edgepreddisthndl,1);
	nedges_traversed=0;

	for (i = 0; i < nlocalverts; ++i)
		for(j = vrowstarts[i];j<vrowstarts[i+1];j++)
			sendedgepreddist(i,j);
	aml_barrier();

	for (i = 0; i < nlocalverts; ++i)
		if(confirmed[i]==0 && pred[i]!=-1)
		printf("Validation Error: path to vertex %ld not confirmed from predecessor %ld with distance %f\n",VERTEX_TO_GLOBAL(my_pe(),i),pred[i],dist[i]),val_errors++;

	aml_long_allsum(&val_errors);
	aml_long_allsum(&nedges_traversed);
	if(nedges_in!=NULL && *nedges_in!=nedges_traversed)
	printf("Validation Error: wrong nedge_traversed %lu (correct number is %lu)\n",*nedges_in,nedges_traversed),val_errors++;
	if(val_errors>0) return 0; else return 1;
}
