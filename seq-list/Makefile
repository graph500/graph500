CC = gcc-4.3 -g -Wall -std=gnu9x -march=native -O3 -DNDEBUG
LDLIBS = -lm -lrt

OBJS = graph500.o graph500-seq-list.o graph500-prng.o graph500-util.o graph500-rmat.o

graph500-seq-list:	$(OBJS)

graph500.o:	graph500.c graph500.h
graph500-seq-list.o:	graph500-seq-list.c graph500.h
graph500-prng.o:	graph500-prng.c graph500-prng.h
graph500-util.o:	graph500-util.c graph500-util.h
graph500-rmat.o:	graph500-rmat.c graph500-rmat.h

.PHONY:	clean
clean:
	rm graph500-seq-list $(OBJS)
