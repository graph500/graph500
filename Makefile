CFLAGS = -g -std=c99
LDLIBS = -lm -lrt

GRAPH500_SOURCES=graph500.c xalloc.c timer.c 
GRAPH500_SEQ_SOURCES=rmat.c verify.c prng.c

BIN=seq-list/seq-list seq-csr/seq-csr

.PHONY: all
all: $(BIN)

seq-list/seq-list: seq-list/seq-list.c $(GRAPH500_SOURCES) $(GRAPH500_SEQ_SOURCES) libgenerator-seq.a
seq-csr/seq-csr: seq-csr/seq-csr.c $(GRAPH500_SOURCES) $(GRAPH500_SEQ_SOURCES) libgenerator-seq.a

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
