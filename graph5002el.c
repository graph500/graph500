/*  Converts from the binary format used in the graph500 to the simpler
 *  edge list format of one <vertex1> <vertex2> pair per line in ASCII
 *  Assumes you can fit the entire graph in memory (the same assumption
 *  the graph500 makes)
 *  Author: Samuel Pollard
 */
#include <stdio.h>
#include <stdlib.h>
#include "xalloc.h" // for xmaloc_large_ext
#include "generator/make_graph.h" // for packed_edge
#include "options.h" // for NBFS_max

static const char* usage =
    "usage: graph5002el <inputedges> <inputroots> <outputedges> <outputroots>";
static struct packed_edge * restrict IJ;

int main (int argc, char **argv)
{
  if (argc != 5) {
    printf ("%s\n", usage);
    return EXIT_FAILURE;
  }

  /* Convert the roots */
  ssize_t sz, sz_read;
  FILE* infp;
  if ((infp  = fopen (argv[2], "r")) == NULL) {
    perror("Can't read file at argv[2]");
    return EXIT_FAILURE;
  }
  fseek (infp, 0L, SEEK_END);
  sz = ftell (infp);
  rewind (infp);
  printf ("Found %ld roots\n", sz / 8); // VERBOSE
  int64_t bfs_root[sz / 8]; // ~64, small enough for the stack
  if ((sz_read = fread (bfs_root, 1, sz, infp)) != sz) {
    fprintf (stderr, "Expected %ld bytes, read %ld\n", sz, sz_read);
    return EXIT_FAILURE;
  }
  fclose (infp);
  FILE* outfp;
  if ((outfp = fopen (argv[4], "w")) == NULL) {
    perror("Can't write file at argv[4]");
    return EXIT_FAILURE;
  }
  int64_t m;
  for (m = 0; m < sz / 8; m++) {
    fprintf (outfp, "%lld\n", bfs_root[m]);
  }
  fclose (outfp);

  /* Convert the edge list */ 
  if ((infp  = fopen (argv[1], "r")) == NULL) {
    perror("Can't read file at argv[1]");
    return EXIT_FAILURE;
  }
  fseek (infp, 0L, SEEK_END);
  sz = ftell (infp);
  rewind (infp);
  /* This is a highly suboptimal way of doing things.
   * If the graph can't fit in RAM, this could be modifed to loop
   * making sure to read in multiples of sizeof(struct packed_edge)
   */
  const size_t chunks = 1;
  IJ = xmalloc_large_ext (sz / chunks);
  printf ("Found %ld edges\n", sz / sizeof(*IJ)); // VERBOSE
  if ((sz_read = fread (IJ, 1, sz, infp)) != sz) {
    fprintf (stderr, "Expected %ld bytes, read %ld\n", sz, sz_read);
    return EXIT_FAILURE;
  }
  fclose (infp);
  if ((outfp = fopen (argv[4], "w")) == NULL) {
    perror("Can't write file at argv[3]");
    return EXIT_FAILURE;
  }
  int64_t i,j;
  for (m = 0; m < sz / sizeof(*IJ); m++) {
    i = get_v0_from_edge(&IJ[m]);
    j = get_v1_from_edge(&IJ[m]);
    fprintf (outfp, "%lld %lld\n", i, j);
  }

  /* Cleanup */
  fclose (outfp);
  xfree_large (IJ);
  return EXIT_SUCCESS;
}
