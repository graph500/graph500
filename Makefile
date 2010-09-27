# -*- Makefile -*-
BUILD_OPENMP = No
BUILD_XMT = No
include make.inc

GRAPH500_SOURCES=graph500.c options.c rmat.c verify.c prng.c xalloc.c timer.c 

BIN=seq-list/seq-list seq-csr/seq-csr

ifeq ($(BUILD_OPENMP), Yes)
BIN += omp-csr/omp-csr
endif

ifeq ($(BUILD_XMT), Yes)
BIN += xmt-csr/xmt-csr
endif

.PHONY: all
all: $(BIN)

seq-list/seq-list: seq-list/seq-list.c $(GRAPH500_SOURCES) #libgenerator-seq.a
seq-csr/seq-csr: seq-csr/seq-csr.c $(GRAPH500_SOURCES) #libgenerator-seq.a

omp-csr/omp-csr: CFLAGS:=$(CFLAGS) $(CFLAGS_OPENMP)
omp-csr/omp-csr: omp-csr/omp-csr.c $(GRAPH500_SOURCES)

xmt-csr/xmt-csr: CFLAGS:=$(CFLAGS) -pl xmt-csr/xmt-csr.pl
#xmt-csr/xmt-csr: CPPFLAGS:=$(CPPFLAGS) -Drestrict=
xmt-csr/xmt-csr: xmt-csr/xmt-csr.c $(GRAPH500_SOURCES)

GENERATOR_OBJS_SEQ=btrd_binomial_distribution.o splittable_mrg.o	\
	mrg_transitions.o graph_generator.o permutation_gen.o		\
	make_graph.o

generator/generator_test_seq: generator/generator_test_seq.c libgenerator-seq.a

libgenerator-seq.a: libgenerator-seq.a($(addprefix generator/,$(GENERATOR_OBJS_SEQ)))
	ranlib libgenerator-seq.a

libgenerator-seq.a($(addprefix generator/,$(GENERATOR_OBJS_SEQ))): CPPFLAGS=-DGRAPH_GENERATOR_SEQ

.PHONY:	clean
clean:
	rm -f libgenerator-seq.a generator/generator_test_seq $(BIN)
