CFLAGS = -fopenmp -std=c99 -g -Wall -Drestrict=__restrict__ -O3 -DGRAPH_GENERATOR_MPI # -g -pg
# CFLAGS = -g -Wall -Drestrict= -DGRAPH_GENERATOR_MPI
LDFLAGS = -fopenmp -g # -g -pg
CPPFLAGS+=-I../support/Random123/include
MPICC = mpicc

all: graph500_mpi_simple graph500_mpi_one_sided graph500_mpi_replicated graph500_mpi_replicated_csc graph500_mpi_custom

GENERATOR_SOURCES = ../generator.c ../globals.c ../prng.c ../xalloc.c
SOURCES = main.c oned_csr.c oned_csc.c utils.c validate.c onesided.c onesided_emul.c
HEADERS = common.h oned_csr.h oned_csc.h redistribute.h mpi_workarounds.h onesided.h

graph500_mpi_simple: bfs_simple.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o graph500_mpi_simple bfs_simple.c $(SOURCES) $(GENERATOR_SOURCES) -lm

graph500_mpi_one_sided: bfs_one_sided.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o graph500_mpi_one_sided bfs_one_sided.c $(SOURCES) $(GENERATOR_SOURCES) -lm

graph500_mpi_replicated: bfs_replicated.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o graph500_mpi_replicated bfs_replicated.c $(SOURCES) $(GENERATOR_SOURCES) -lm

graph500_mpi_replicated_csc: bfs_replicated_csc.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o graph500_mpi_replicated_csc bfs_replicated_csc.c $(SOURCES) $(GENERATOR_SOURCES) -lm

graph500_mpi_custom: bfs_custom.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o graph500_mpi_custom bfs_custom.c $(SOURCES) $(GENERATOR_SOURCES) -lm

clean:
	-rm -f graph500_mpi_* *.o *.a
