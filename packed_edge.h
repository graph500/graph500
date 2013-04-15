#if !defined(PACKED_EDGE_HEADER_)
#define PACKED_EDGE_HEADER_

#ifdef GENERATOR_USE_PACKED_EDGE_TYPE

typedef struct packed_edge {
  uint8_t data[5+5+1];
} packed_edge;

static inline int64_t get_v0_from_edge(const packed_edge* p) {
  int64_t out = 0;
  const uint8_t * restrict data = p->data;
  for (int k = 0; k < 5; ++k)
    out |= (int64_t)(data[k] << (8*k));
  return out;
}

static inline int64_t get_v1_from_edge(const packed_edge* p) {
  int64_t out = 0;
  const uint8_t * restrict data = p->data;
  for (int k = 0; k < 5; ++k)
    out |= (int64_t)(data[5+k] << (8*k));
  return out;
}

static inline uint8_t get_w_from_edge(const packed_edge* p) {
  return p->data[11];
}

static inline void write_edge(packed_edge* p, int64_t v0, int64_t v1,
			      uint8_t w)
{
  uint8_t * restrict data = p->data;
  uint64_t uv0 = v0;
  uint64_t uv1 = v1;
  for (int k = 0; k < 5; ++k) {
    data[k] = (uint8_t)(uv0 & (uint64_t)0xFF);
    data[5+k] = (uint8_t)(uv1 & (uint64_t)0xFF);
    uv0 >>= 8;
    uv1 >>= 8;
  }
  data[11] = w;
}

#else

typedef struct packed_edge {
  int64_t v0;
  int64_t v1;
  uint8_t w;
} packed_edge;

static inline int64_t get_v0_from_edge(const packed_edge* p) {
  return p->v0;
}

static inline int64_t get_v1_from_edge(const packed_edge* p) {
  return p->v1;
}

static inline uint8_t get_w_from_edge(const packed_edge* p) {
  return p->w;
}

static inline void write_edge(packed_edge* p, int64_t v0, int64_t v1,
			      uint8_t w)
{
  p->v0 = v0;
  p->v1 = v1;
  p->w = w;
}

#endif

#endif /* PACKED_EDGE_HEADER_ */
