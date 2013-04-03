#if !defined(GENERATOR_HEADER_)
#define GENERATOR_HEADER_

void make_edge (int64_t, int64_t * restrict, int64_t * restrict,
		uint8_t * restrict);
void edge_list (int64_t * restrict, int64_t * restrict, uint8_t * restrict,
		const int64_t, const int64_t);

/* To go away... */
void make_graph (packed_edge *);
void packed_edge_list (packed_edge *, const int64_t, const int64_t);

#endif /* GENERATOR_HEADER_ */
