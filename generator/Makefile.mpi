CC = mpicc
CFLAGS = -std=c99 -O3 -DGRAPH_GENERATOR_MPI -DGRAPHGEN_DISTRIBUTED_MEMORY -DNDEBUG # -g -pg
# CFLAGS = -std=c99 -DGRAPH_GENERATOR_MPI -DGRAPHGEN_DISTRIBUTED_MEMORY -g
LDFLAGS = -O3
# LDFLAGS = -g
MPICC = mpicc

GENERATOR_SOURCES = graph_generator.c make_graph.c splittable_mrg.c utils.c
GENERATOR_HEADERS = graph_generator.h make_graph.h mod_arith_32bit.h mod_arith_64bit.h mod_arith.h mod_arith_xmt.h splittable_mrg.h utils.h user_settings.h mrg_transitions.c

all: generator_test_mpi

generator_test_mpi: generator_test_mpi.c $(GENERATOR_SOURCES) $(GENERATOR_HEADERS)
	$(MPICC) $(CFLAGS) $(LDFLAGS) -o generator_test_mpi generator_test_mpi.c $(GENERATOR_SOURCES) -lm

clean:
	-rm -f generator_test_mpi
