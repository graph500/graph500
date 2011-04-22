#if !defined(KRONECKER_HEADER_)
#define KRONECKER_HEADER_

#include "generator/graph_generator.h"

void kronecker_edgelist (struct packed_edge *IJ, int64_t nedge, int64_t SCALE,
			 double A, double B, double C);

#endif /* KRONECKER_HEADER_ */
