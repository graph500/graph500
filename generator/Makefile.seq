AR = ar
RANLIB = ranlib
CC = cc
CFLAGS = -g -Wall -Drestrict=__restrict__ -O3 -DNDEBUG -ffast-math -DGRAPH_GENERATOR_SEQ # -g -pg
# CFLAGS = -g -Wall -Drestrict=__restrict__
LDFLAGS = -g # -g -pg

all: libgraph_generator_seq.a generator_test_seq

generator_test_seq: generator_test_seq.c libgraph_generator_seq.a
	$(CC) $(CFLAGS) $(LDFLAGS) -o generator_test_seq generator_test_seq.c -L. -lgraph_generator_seq -lm

libgraph_generator_seq.a: btrd_binomial_distribution.o splittable_mrg.o mrg_transitions.o graph_generator.o permutation_gen.o make_graph.o utils.o scramble_edges.o
	$(AR) cruv libgraph_generator_seq.a btrd_binomial_distribution.o splittable_mrg.o mrg_transitions.o graph_generator.o permutation_gen.o make_graph.o utils.o scramble_edges.o
	$(RANLIB) libgraph_generator_seq.a

clean:
	-rm -f generator_test_seq libgraph_generator_seq.a *.o
