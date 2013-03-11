# -*- Makefile -*-
# Copyright 2010-2011,  Georgia Institute of Technology, USA.
# See COPYING for license.
BUILD_OPENMP = No
BUILD_XMT = No
include make.inc

GRAPH500_SOURCES=graph500.c options.c rmat.c kronecker.c verify.c prng.c \
	xalloc.c timer.c 

MAKE_EDGELIST_SOURCES=make-edgelist.c options.c rmat.c kronecker.c prng.c \
	xalloc.c timer.c 

BIN=seq-list/seq-list seq-csr/seq-csr make-edgelist

ifeq ($(BUILD_OPENMP), Yes)
BIN += omp-csr/omp-csr
endif

ifeq ($(BUILD_MPI), Yes)
BIN += mpi/graph500_mpi_simple
endif

ifeq ($(BUILD_XMT), Yes)
BIN = xmt-csr/xmt-csr xmt-csr-local/xmt-csr-local
endif

GENERATOR_SRCS=splittable_mrg.c	\
	graph_generator.c \
	make_graph.c utils.c

.PHONY: all
all: $(BIN)

CPPFLAGS += -I./generator

make-edgelist: CFLAGS:=$(CFLAGS) $(CFLAGS_OPENMP)
make-edgelist:	make-edgelist.c options.c rmat.c kronecker.c prng.c \
	xalloc.c timer.c $(addprefix generator/,$(GENERATOR_SRCS))

seq-list/seq-list: seq-list/seq-list.c $(GRAPH500_SOURCES) \
	$(addprefix generator/,$(GENERATOR_SRCS))
seq-csr/seq-csr: seq-csr/seq-csr.c $(GRAPH500_SOURCES) \
	$(addprefix generator/,$(GENERATOR_SRCS))

omp-csr/omp-csr-old: CFLAGS:=$(CFLAGS) $(CFLAGS_OPENMP)
omp-csr/omp-csr-old: omp-csr/omp-csr-old.c $(GRAPH500_SOURCES) \
	$(addprefix generator/,$(GENERATOR_SRCS))

omp-csr/omp-csr: CFLAGS:=$(CFLAGS) $(CFLAGS_OPENMP)
omp-csr/omp-csr: omp-csr/omp-csr.c omp-csr/bitmap.h $(GRAPH500_SOURCES) \
	$(addprefix generator/,$(GENERATOR_SRCS))

xmt-csr/xmt-csr: CFLAGS:=$(CFLAGS) -pl xmt-csr/xmt-csr.pl
xmt-csr/xmt-csr: xmt-csr/xmt-csr.c $(GRAPH500_SOURCES) \
	$(addprefix generator/,$(GENERATOR_SRCS))

xmt-csr-local/xmt-csr-local: CFLAGS:=$(CFLAGS) -pl xmt-csr-local/xmt-csr-local.pl
xmt-csr-local/xmt-csr-local: xmt-csr-local/xmt-csr-local.c $(GRAPH500_SOURCES) \
	$(addprefix generator/,$(GENERATOR_SRCS))

generator/generator_test_seq: generator/generator_test_seq.c $(GENERATOR_SRCS)

generator/generator_test_omp: generator/generator_test_omp.c $(GENERATOR_SRCS)

mpi/graph500_mpi_simple mpi/graph500_mpi_one_sided mpi/graph500_mpi_replicated:
	$(MAKE) -C mpi

.PHONY:	clean
clean:
	rm -f generator/generator_test_seq generator/generator_test_omp \
		$(BIN)
	-$(MAKE) -C mpi clean
