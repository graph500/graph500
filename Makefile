# -*- Makefile -*-
# Copyright 2010,  Georgia Institute of Technology, USA.
# See COPYING for license.
BUILD_OPENMP = No
BUILD_XMT = No
include make.inc

GRAPH500_SOURCES=graph500.c options.c rmat.c kronecker.c verify.c prng.c \
	xalloc.c timer.c 

BIN=seq-list/seq-list seq-csr/seq-csr

ifeq ($(BUILD_OPENMP), Yes)
BIN += omp-csr/omp-csr
endif

ifeq ($(BUILD_MPI), Yes)
BIN += mpi/graph500_mpi_simple
endif

ifeq ($(BUILD_XMT), Yes)
BIN = xmt-csr/xmt-csr xmt-csr-local/xmt-csr-local
endif

GENERATOR_OBJS_SEQ=btrd_binomial_distribution.o splittable_mrg.o	\
	mrg_transitions.o graph_generator.o permutation_gen.o		\
	make_graph.o scramble_edges.o utils.o

.PHONY: all
all: $(BIN)

seq-list/seq-list: seq-list/seq-list.c $(GRAPH500_SOURCES) libgenerator-seq.a
seq-csr/seq-csr: seq-csr/seq-csr.c $(GRAPH500_SOURCES) libgenerator-seq.a

omp-csr/omp-csr: CFLAGS:=$(CFLAGS) $(CFLAGS_OPENMP)
omp-csr/omp-csr: omp-csr/omp-csr.c $(GRAPH500_SOURCES) libgenerator-omp.a

xmt-csr/xmt-csr: CFLAGS:=$(CFLAGS) -pl xmt-csr/xmt-csr.pl
xmt-csr/xmt-csr: xmt-csr/xmt-csr.c $(GRAPH500_SOURCES) \
	$(addprefix generator/,$(patsubst %.o,%.c,$(GENERATOR_OBJS_SEQ)))

xmt-csr-local/xmt-csr-local: CFLAGS:=$(CFLAGS) -pl xmt-csr-local/xmt-csr-local.pl
xmt-csr-local/xmt-csr-local: xmt-csr-local/xmt-csr-local.c $(GRAPH500_SOURCES) \
	$(addprefix generator/,$(patsubst %.o,%.c,$(GENERATOR_OBJS_SEQ)))

generator/generator_test_seq: generator/generator_test_seq.c libgenerator-seq.a

generator/generator_test_omp: generator/generator_test_omp.c libgenerator-omp.a

libgenerator-seq.a: libgenerator-seq.a($(addprefix generator/,$(GENERATOR_OBJS_SEQ)))
	ranlib libgenerator-seq.a

libgenerator-seq.a($(addprefix generator/,$(GENERATOR_OBJS_SEQ))): CFLAGS:=$(CFLAGS) $(CFLAGS_OPENMP)
libgenerator-seq.a($(addprefix generator/,$(GENERATOR_OBJS_SEQ))): CPPFLAGS=-DGRAPH_GENERATOR_SEQ

libgenerator-omp.a: libgenerator-omp.a($(addprefix generator/,$(GENERATOR_OBJS_SEQ)))
	ranlib libgenerator-omp.a

libgenerator-omp.a($(addprefix generator/,$(GENERATOR_OBJS_SEQ))): CPPFLAGS=-DGRAPH_GENERATOR_OMP

mpi/graph500_mpi_simple mpi/graph500_mpi_one_sided:
	$(MAKE) -C mpi

.PHONY:	clean
clean:
	rm -f libgenerator-omp.a libgenerator-seq.a \
		generator/generator_test_seq generator/generator_test_omp \
		$(BIN)
	-$(MAKE) -C mpi clean
