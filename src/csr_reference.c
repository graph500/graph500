/* Copyright (c) 2011-2017 Graph500 Steering Committee
   All rights reserved.
   Developed by:                Anton Korzh anton@korzh.us
                                Graph500 Steering Committee
                                http://www.graph500.org
   New code under University of Illinois/NCSA Open Source License
   see license.txt or https://opensource.org/licenses/NCSA
*/

// Graph500: Kernel 1: CRS construction
// Simple two-pass CRS construction using Active Messages

#include "common.h"
#include "csr_reference.h"
#include "aml.h"
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <search.h>

int64_t nverts_known = 0;
int *degrees;
int64_t *column;
float *weights;
extern oned_csr_graph g; //from bfs_reference for isisolated function

//this function is needed for roots generation
int isisolated(int64_t v) {
	if(my_pe()==VERTEX_OWNER(v)) return (g.rowstarts[VERTEX_LOCAL(v)]==g.rowstarts[VERTEX_LOCAL(v)+1]);
	return 0; //locally no evidence, allreduce required
}

void halfedgehndl(int from,void* data,int sz)
{  degrees[*(int*)data]++; }

void fulledgehndl(int frompe,void* data,int sz) {
	int vloc = *(int*)data;
	int64_t gtgt = *((int64_t*)(data+4));
	SETCOLUMN(degrees[vloc]++,gtgt);
#ifdef SSSP
	float w = ((float*)data)[3];
	weights[degrees[vloc]-1]=w;
#endif
}

void send_half_edge (int64_t src,int64_t tgt) {
	int pe=VERTEX_OWNER(src);
	int vloc=VERTEX_LOCAL(src);
	aml_send(&vloc,1,4,pe);
	if(tgt>=nverts_known) nverts_known=tgt+1;
}
#ifdef SSSP
void send_full_edge (int64_t src,int64_t tgt,float w) {
	int pe=VERTEX_OWNER(src);
	int vloc[4];
	vloc[0]=VERTEX_LOCAL(src);
	memcpy(vloc+1,&tgt,8);
	memcpy(vloc+3,&w,4);
	aml_send(vloc,1,16,pe);
}
#else
void send_full_edge (int64_t src,int64_t tgt) {
	int pe=VERTEX_OWNER(src);
	int vloc[3];
	vloc[0]=VERTEX_LOCAL(src);
	memcpy(vloc+1,&tgt,8);
	aml_send(vloc,1,12,pe);
}
#endif

void convert_graph_to_oned_csr(const tuple_graph* const tg, oned_csr_graph* const g) {
	g->tg = tg;

	size_t i;

	int64_t nvert=tg->nglobaledges/2;
	nvert/=num_pes();
	nvert+=1;
	degrees=xcalloc(nvert,sizeof(int));

	aml_register_handler(halfedgehndl,1);
	// First pass : calculate degrees of each vertex
	ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize,wbuf) {
		ptrdiff_t j;
		for (j = 0; j < bufsize; ++j) {
			int64_t v0 = get_v0_from_edge(&buf[j]);
			int64_t v1 = get_v1_from_edge(&buf[j]);
			if(v0==v1) continue;
			send_half_edge(v0, v1);
			send_half_edge(v1, v0);
		}
		aml_barrier();
	} ITERATE_TUPLE_GRAPH_END;

	int64_t nglobalverts = 0;
	aml_long_allmax(&nverts_known);
	nglobalverts=nverts_known+1;
	g->nglobalverts = nglobalverts;
	size_t nlocalverts = VERTEX_LOCAL(nglobalverts + num_pes() - 1 - my_pe());
	g->nlocalverts = nlocalverts;

	//graph stats printing
#ifdef DEBUGSTATS
	long maxdeg=0,isolated=0,totaledges=0,originaledges;
	long maxlocaledges,minlocaledges;
	for(i=0;i<g->nlocalverts;i++) {
		long deg = degrees[i];
		totaledges+=deg;
		if(maxdeg<deg) maxdeg=deg;
		if(!deg) isolated++;
	}
	originaledges=totaledges;
	maxlocaledges=totaledges;
	minlocaledges=totaledges;
	aml_long_allmax(&maxdeg);
	aml_long_allsum(&isolated);
	aml_long_allsum(&totaledges);
	aml_long_allmin(&minlocaledges);
	aml_long_allmax(&maxlocaledges);
	long averageedges = totaledges/num_pes();
	double disbalance = (double)(maxlocaledges-minlocaledges)/(double)averageedges * 100.0;
	if(!my_pe()) printf("\n maxdeg %lld verts %lld, isolated %lld edges %lld\n\t A max %ld min %ld ave %ld delta %ld percent %3.2f\n ",
			maxdeg,g->nglobalverts,isolated,totaledges,maxlocaledges,minlocaledges,averageedges,maxlocaledges-minlocaledges,disbalance);

	// finished stats printing

	g->notisolated=g->nglobalverts-isolated;
#endif
	unsigned int *rowstarts = xmalloc((nlocalverts + 1) * sizeof(int));
	g->rowstarts = rowstarts;

	rowstarts[0] = 0;
	for (i = 0; i < nlocalverts; ++i) {
		rowstarts[i + 1] = rowstarts[i] + (i >= nlocalverts ? 0 : degrees[i]);
		degrees[i] = rowstarts[i];
	}

	size_t nlocaledges = rowstarts[nlocalverts];
	g->nlocaledges = nlocaledges;

	int64_t colalloc = BYTES_PER_VERTEX*nlocaledges;
	colalloc += (4095);
	colalloc /= 4096;
	colalloc *= 4096;
	column = xmalloc(colalloc);
	aml_barrier();
#ifdef SSSP
	weights = xmalloc(4*nlocaledges);
	g->weights = weights;
	aml_barrier();
#endif
	//long allocatededges=colalloc;
	g->column = column;

	aml_register_handler(fulledgehndl,1);
	//Next pass , actual data transfer: placing edges to its places in column and hcolumn
	ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize,wbuf) {
		ptrdiff_t j;
		for (j = 0; j < bufsize; ++j) {
			int64_t v0 = get_v0_from_edge(&buf[j]);
			int64_t v1 = get_v1_from_edge(&buf[j]);
			if(v0==v1) continue;
#ifdef SSSP
			send_full_edge(v0, v1,wbuf[j]);
			send_full_edge(v1, v0,wbuf[j]);
#else				
			send_full_edge(v0, v1);
			send_full_edge(v1, v0);
#endif
		}
		aml_barrier();
	} ITERATE_TUPLE_GRAPH_END;

	free(degrees);
}

void free_oned_csr_graph(oned_csr_graph* const g) {
	if (g->rowstarts != NULL) {free(g->rowstarts); g->rowstarts = NULL;}
	if (g->column != NULL) {free(g->column); g->column = NULL;}
#ifdef SSSP
	if (g->weights != NULL) {free(g->weights); g->weights = NULL;}
#endif
}
