Graph500 Benchmark: MPI Reference Implementation (BFS and SSSP kernels) 
Version 3.0.0rc1 
Jeremiah Willcock and Andrew Lumsdaine, Anton Korzh

This directory contains an MPI-based reference implementation for the Graph500
benchmark.  The Makefile will automatically build the supporting code, along
with two reference implementations one doing both breadth-first search and
single source shortest path, other just doing BFS.

Makefile may need to be modified with your MPI location and compiler flags.
The generated programs, named graph500_mpi_*, take up to two parameters: the
scale of the problem and the edge factor.  The problem scale is the logarithm,
base 2, of the number of vertices in the graph; only graphs with power-of-2
vertex counts are supported without source code modification.  The edge factor
is the ratio of the number of edges to the number of vertices; i.e., it is half
the average vertex degree in the graph.  The scale parameter is mandatory; the
edge factor is optional and defaults to 16 (the value specified by the
benchmark).  Running any of the graph500_* programs without any arguments will
produce a usage message.

The raw graph tuples produced by the graph generator can be stored either in
memory or in a file (using MPI file I/O).  If the TMPFILE environment variable
is set, it should point to a filename (in the sense of section 13.2.1 of the
MPI 2.2 standard) whose corresponding file does not exist and which points to a
filesystem that is globally accessible and consistent on all ranks and has
enough storage space for the graph data.  The amount of space required is given
in <URL:http://www.graph500.org/Specifications.html#tbl:classes> or can be
computed as 256*(2**SCALE) bytes.  If GENERATOR_USE_PACKED_EDGE_TYPE is
#define'd in ../generator/user_settings.h, the space required is reduced to
192*(2**SCALE) bytes. The bfs-sssp binary will also produce a separate weights
file whose name would be TMPFILE appended by ".weights".  Edges file is exactly
the same as it was for reference 2.1.4.  In case edges file is present but no
weights file exists both would be recreated.

The code is written in C; the code compiles with GCC's default gnu89 language
setting, but should be valid C99 and C++ (except for the use of a few C99
headers).  The main non-C89 features used are variable declarations after
statements in a block and the <stdint.h> and <inttypes.h> headers.  The code
assumes that your MPI implementation is compliant with version 3.0 of the MPI
standard; it uses some MPI non-blocking collectives not present in MPI 2.*. All
communication is implemented through AML (Active Messages Library) which
implements asynchronous active messages and implements transparent message
coalescing (or aggregation) and software two-step routing for efficient
performance of modern multicore nodes. No need of using OpenMP or Hybrid mode.


A template for writing your own BFS using the reference data structures and
infrastructure is in bfs_custom.c.  You can either modify that file in place or
copy it (adjusting the Makefile) to create your own version.  The documentation
for what data structures are available and how to use them is in comments in
bfs_custom.c.  The reference implementation also contains code to convert from
a distributed list of graph edges into a distributed compressed sparse row data
structure, as well as code for timing the BFS run, validating the correctness
of the results, and printing the timings in the Graph500-required format.


------------------------------------------------------------------------ 
Flags and variables:

- compiling without SSSP macro defined build only BFS code
- compiling with SSSP macro defined build code which runs BFS and then SSSP
- for latter its possible do skip BFS by having env variable SKIP_BFS
- Validation can be skipped by having env variable SKIP_VALIDATION, but given
  new validation is much faster,less memory consuming (its only 30% slower then
BFS run and faster then SSSP run ) its recommended to have validation always
enabled, unless you really need all your memory for BFS/SSSP run
- macro REUSE_CSR_FOR_VALIDATION if enabled would reuse reference CSR generated
  by kernel 1 convert for validation thus greatly saving memory (validation
would require almost no additional memory to proceed, otherwise it builds it's
own reference CRS)
- macro DEBUGSTATS when enabled gives nice information about the traversed
  graph levels

Troubleshooting:

-If number of processes per node is not a power of 2 you would have to add
macro definition -DPROCS_PER_NODE_NOT_POWER_OF_TWO to CFLAGS in a src/Makefile. 
Macro SIZE_MUST_BE_A_POWER_OF_TWO would be undefined automatically. Experiments
normally give better performance if you use power of twos: using 16 cores on 24-core
node is always faster.

-If you are running on non-power of two of nodes, you would need to comment out
definition of a macro SIZE_MUST_BE_A_POWER_OF_TWO in common.h. But beware it's
known that running on a non power of two nodes slows you down by ~20%

- AML is using decent amount of memory per node for message coalescing.  If
  thats a problem for your machine you should reduce AGGR macro in aml/aml.c
Memory consumed by those buffers is NNODES*NPESPERNODE*AGGR on each node, so
its NPES*AGGR for whatever NPESPERNODE you having.  .So if you are running 8192
processes it would allocate currently 256 MB per node, for 32768 processes it
would be 1GB per node. For 128k processes total it would be 4GB per node.  So
you dont need to touch that parameter unless you have a really big machine,
which hopefully has nice custom interconnect and smaller buffers would work
fine for you.  Default paramters work well both on slow Ethernets and old and
new Infinibands.
------------------------------------------------------------------------
