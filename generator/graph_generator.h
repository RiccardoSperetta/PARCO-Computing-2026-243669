/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#ifndef GRAPH_GENERATOR_H
#define GRAPH_GENERATOR_H

#include "user_settings.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef GENERATOR_USE_PACKED_EDGE_TYPE
#define GENERATOR_USE_PACKED_EDGE_TYPE 0
#endif

#if GENERATOR_USE_PACKED_EDGE_TYPE 0

/**
 * @struct packed_edge
 * @brief A structure representing a packed edge in a graph.
 *
 * This structure is designed to efficiently store edges in a graph by packing
 * vertex identifiers into a compact format. It uses a combination of low and high
 * bits to represent two vertices of an edge.
 *
 * @var packed_edge::v0_low
 * The lower 32 bits of the first vertex identifier (v0). This allows for a 
 * representation of vertex indices that fit within the limits of a 32-bit 
 * unsigned integer.
 *
 * @var packed_edge::v1_low
 * The lower 32 bits of the second vertex identifier (v1). Similar to v0_low, 
 * this allows for efficient storage of the second vertex's index.
 *
 * @var packed_edge::high
 * The higher bits of the edge representation, which can be used to store 
 * additional information about the edge, such as the higher bits of v1. 
 * This allows for a larger range of vertex identifiers to be represented 
 * within a single packed_edge structure.
 *
 * The use of both v0_low and v1_low allows for the representation of two 
 * vertices in a single structure while maintaining the ability to access 
 * each vertex's identifier efficiently.
 */
typedef struct packed_edge {
  uint32_t v0_low;
  uint32_t v1_low;
  uint32_t high; /* v1 in high half, v0 in low half */
} packed_edge;

static inline int64_t get_v0_from_edge(const packed_edge* p) {
  return (p->v0_low | ((int64_t)((int16_t)(p->high & 0xFFFF)) << 32));
}

static inline int64_t get_v1_from_edge(const packed_edge* p) {
  return (p->v1_low | ((int64_t)((int16_t)(p->high >> 16)) << 32));
}

static inline void write_edge(packed_edge* p, int64_t v0, int64_t v1) {
  p->v0_low = (uint32_t)v0;
  p->v1_low = (uint32_t)v1;
  p->high = (uint32_t)(((v0 >> 32) & 0xFFFF) | (((v1 >> 32) & 0xFFFF) << 16));
}

#else

typedef struct packed_edge {
  int64_t v0;
  int64_t v1;
} packed_edge;

static inline int64_t get_v0_from_edge(const packed_edge* p) {
  return p->v0;
}

static inline int64_t get_v1_from_edge(const packed_edge* p) {
  return p->v1;
}

static inline void write_edge(packed_edge* p, int64_t v0, int64_t v1) {
  p->v0 = v0;
  p->v1 = v1;
}

#endif

/* Generate a range of edges (from start_edge to end_edge of the total graph),
 * writing into elements [0, end_edge - start_edge) of the edges array.  This
 * code is parallel on OpenMP and XMT; it must be used with
 * separately-implemented SPMD parallelism for MPI. */
void generate_kronecker_range(
       const uint_fast32_t seed[5] /* All values in [0, 2^31 - 1) */,
       int logN /* In base 2 */,
       int64_t start_edge, int64_t end_edge /* Indices (in [0, M)) for the edges to generate */,
       packed_edge* edges /* Size >= end_edge - start_edge */
#ifdef SSSP
       ,float* weights
#endif
);

#ifdef __cplusplus
}
#endif

#endif /* GRAPH_GENERATOR_H */
