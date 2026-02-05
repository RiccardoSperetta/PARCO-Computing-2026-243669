/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include <math.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "user_settings.h"
#include "splittable_mrg.h"
#include "graph_generator.h"

/**
 * @file generator.c
 * @brief This file contains the implementation of the generator function.
 *
 * The generator function is designed to produce a sequence of values based on
 * the specified parameters. It can be used in various applications where
 * a series of generated values is required.
 *
 * @note To compile this generator for use, run the following command in your terminal:
 *       gcc -o generator generator.c
 *
 * @usage Example usage of the generator function can be found in the main function
 *        or in the accompanying documentation.
 */

/* Initiator settings: for faster random number generation, the initiator
 * probabilities are defined as fractions (a = INITIATOR_A_NUMERATOR /
 * INITIATOR_DENOMINATOR, b = c = INITIATOR_BC_NUMERATOR /
 * INITIATOR_DENOMINATOR, d = 1 - a - b - c. */
#define INITIATOR_A_NUMERATOR 5700
#define INITIATOR_BC_NUMERATOR 1900
#define INITIATOR_DENOMINATOR 10000

/* If this macro is defined to a non-zero value, use SPK_NOISE_LEVEL /
 * INITIATOR_DENOMINATOR as the noise parameter to use in introducing noise
 * into the graph parameters.  The approach used is from "A Hitchhiker's Guide
 * to Choosing Parameters of Stochastic Kronecker Graphs" by C. Seshadhri, Ali
 * Pinar, and Tamara G. Kolda (http://arxiv.org/abs/1102.5046v1), except that
 * the adjustment here is chosen based on the current level being processed
 * rather than being chosen randomly. */
#define SPK_NOISE_LEVEL 0
/* #define SPK_NOISE_LEVEL 1000 -- in INITIATOR_DENOMINATOR units */

static int generate_4way_bernoulli(mrg_state* st, int level, int nlevels) {
#if SPK_NOISE_LEVEL == 0
  /* Avoid warnings */
  (void)level;
  (void)nlevels;
#endif
  /* Generate a pseudorandom number in the range [0, INITIATOR_DENOMINATOR)
   * without modulo bias. */
  static const uint32_t limit = (UINT32_C(0x7FFFFFFF) % INITIATOR_DENOMINATOR);
  uint32_t val = mrg_get_uint_orig(st);
  if (/* Unlikely */ val < limit) {
    do {
      val = mrg_get_uint_orig(st);
    } while (val < limit);
  }
#if SPK_NOISE_LEVEL == 0
  int spk_noise_factor = 0;
#else
  int spk_noise_factor = 2 * SPK_NOISE_LEVEL * level / nlevels - SPK_NOISE_LEVEL;
#endif
  unsigned int adjusted_bc_numerator = (unsigned int)(INITIATOR_BC_NUMERATOR + spk_noise_factor);
  val %= INITIATOR_DENOMINATOR;
  if (val < adjusted_bc_numerator) return 1;
  val = (uint32_t)(val - adjusted_bc_numerator);
  if (val < adjusted_bc_numerator) return 2;
  val = (uint32_t)(val - adjusted_bc_numerator);
#if SPK_NOISE_LEVEL == 0
  if (val < INITIATOR_A_NUMERATOR) return 0;
#else
  if (val < INITIATOR_A_NUMERATOR * (INITIATOR_DENOMINATOR - 2 * INITIATOR_BC_NUMERATOR) / (INITIATOR_DENOMINATOR - 2 * adjusted_bc_numerator)) return 0;
#endif
#if SPK_NOISE_LEVEL == 0
  /* Avoid warnings */
  (void)level;
  (void)nlevels;
#endif
  return 3;
}

/* Reverse bits in a number; this should be optimized for performance
 * (including using bit- or byte-reverse intrinsics if your platform has them).
 * */
static inline uint64_t bitreverse(uint64_t x) {
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3)
#define USE_GCC_BYTESWAP /* __builtin_bswap* are in 4.3 but not 4.2 */
#endif

#ifdef FAST_64BIT_ARITHMETIC

  /* 64-bit code */
#ifdef USE_GCC_BYTESWAP
  x = __builtin_bswap64(x);
#else
  x = (x >> 32) | (x << 32);
  x = ((x >> 16) & UINT64_C(0x0000FFFF0000FFFF)) | ((x & UINT64_C(0x0000FFFF0000FFFF)) << 16);
  x = ((x >>  8) & UINT64_C(0x00FF00FF00FF00FF)) | ((x & UINT64_C(0x00FF00FF00FF00FF)) <<  8);
#endif
  x = ((x >>  4) & UINT64_C(0x0F0F0F0F0F0F0F0F)) | ((x & UINT64_C(0x0F0F0F0F0F0F0F0F)) <<  4);
  x = ((x >>  2) & UINT64_C(0x3333333333333333)) | ((x & UINT64_C(0x3333333333333333)) <<  2);
  x = ((x >>  1) & UINT64_C(0x5555555555555555)) | ((x & UINT64_C(0x5555555555555555)) <<  1);
  return x;

#else

  /* 32-bit code */
  uint32_t h = (uint32_t)(x >> 32);
  uint32_t l = (uint32_t)(x & UINT32_MAX);
#ifdef USE_GCC_BYTESWAP
  h = __builtin_bswap32(h);
  l = __builtin_bswap32(l);
#else
  h = (h >> 16) | (h << 16);
  l = (l >> 16) | (l << 16);
  h = ((h >> 8) & UINT32_C(0x00FF00FF)) | ((h & UINT32_C(0x00FF00FF)) << 8);
  l = ((l >> 8) & UINT32_C(0x00FF00FF)) | ((l & UINT32_C(0x00FF00FF)) << 8);
#endif
  h = ((h >> 4) & UINT32_C(0x0F0F0F0F)) | ((h & UINT32_C(0x0F0F0F0F)) << 4);
  l = ((l >> 4) & UINT32_C(0x0F0F0F0F)) | ((l & UINT32_C(0x0F0F0F0F)) << 4);
  h = ((h >> 2) & UINT32_C(0x33333333)) | ((h & UINT32_C(0x33333333)) << 2);
  l = ((l >> 2) & UINT32_C(0x33333333)) | ((l & UINT32_C(0x33333333)) << 2);
  h = ((h >> 1) & UINT32_C(0x55555555)) | ((h & UINT32_C(0x55555555)) << 1);
  l = ((l >> 1) & UINT32_C(0x55555555)) | ((l & UINT32_C(0x55555555)) << 1);
  return ((uint64_t)l << 32) | h; /* Swap halves */

#endif
}

/* Apply a permutation to scramble vertex numbers; a randomly generated
 * permutation is not used because applying it at scale is too expensive. */
static inline int64_t scramble(int64_t v0, int lgN, uint64_t val0, uint64_t val1) {
  uint64_t v = (uint64_t)v0;
  v += val0 + val1;
  v *= (val0 | UINT64_C(0x4519840211493211));
  v = (bitreverse(v) >> (64 - lgN));
  assert ((v >> lgN) == 0);
  v *= (val1 | UINT64_C(0x3050852102C843A5));
  v = (bitreverse(v) >> (64 - lgN));
  assert ((v >> lgN) == 0);
  return (int64_t)v;
}

/* Make a single graph edge using a pre-set MRG state. */
static
void make_one_edge(int64_t nverts, int level, int lgN, mrg_state* st, packed_edge* result, uint64_t val0, uint64_t val1) {
  int64_t base_src = 0, base_tgt = 0;
  while (nverts > 1) {
    int square = generate_4way_bernoulli(st, level, lgN);
    int src_offset = square / 2;
    int tgt_offset = square % 2;
    assert (base_src <= base_tgt);
    if (base_src == base_tgt) {
      /* Clip-and-flip for undirected graph */
      if (src_offset > tgt_offset) {
        int temp = src_offset;
        src_offset = tgt_offset;
        tgt_offset = temp;
      }
    }
    nverts /= 2;
    ++level;
    base_src += nverts * src_offset;
    base_tgt += nverts * tgt_offset;
  }
  write_edge(result,
             scramble(base_src, lgN, val0, val1),
             scramble(base_tgt, lgN, val0, val1));
}


/* ====================================
actual generator function
==================================== */
/* Generate a range of edges (from start_edge to end_edge of the total graph),
 * writing into elements [0, end_edge - start_edge) of the edges array.  This
 * code is parallel on OpenMP and XMT; it must be used with
 * separately-implemented SPMD parallelism for MPI. */
/**
 * @brief Generates a range of edges for a Kronecker graph using the MRG random number generator.
 * 
 * This function generates edges between start_edge and end_edge for a Kronecker graph with
 * 2^logN vertices. The function is parallelizable and can optionally generate weights for
 * each edge when SSSP is defined at compile time.
 * 
 * @param seed Array of 5 seed values for the MRG random number generator.
 *             Each value must be in range [0, 2^31 - 1) and not all zero.
 * @param logN Logarithm base 2 of the number of vertices (vertices = 2^logN).
 * @param start_edge Starting edge index (inclusive) to generate.
 * @param end_edge Ending edge index (exclusive) to generate.
 * @param edges Output array to store generated edges. Must have at least (end_edge - start_edge) capacity.
 * @param weights Optional output array for edge weights (only if compiled with -DSSSP).
 *                Must have at least (end_edge - start_edge) capacity.
 * 
 * @note This function is OpenMP parallelizable when compiled with -fopenmp.
 * @note Weights are only generated if SSSP preprocessor flag is defined at compile time.
 * 
 * @example
 * // Generate edges 1000 to 2000 for a graph with 2^20 vertices
 * uint_fast32_t seed[5] = {12345, 67890, 11111, 22222, 33333};
 * packed_edge edges[1000];
 * float weights[1000];  // Only needed if compiled with -DSSSP
 * 
 * generate_kronecker_range(seed, 20, 1000, 2000, edges, weights);
 * // Result: edges[] contains 1000 randomly generated Kronecker edges
 * //         weights[] contains 1000 random float weights in range [0, 1)
 */
void generate_kronecker_range(
       const uint_fast32_t seed[5] /* All values in [0, 2^31 - 1), not all zero */,
       int logN /* In base 2 */,
       int64_t start_edge, int64_t end_edge,
       packed_edge* edges
       ) {
  mrg_state state;
  int64_t nverts = (int64_t)1 << logN;
  int64_t ei;

  /* Initialize the base RNG state with the provided seed */
  mrg_seed(&state, seed);

  /* Generate two 64-bit scrambling values by advancing the RNG state.
   * These values are used to permute vertex indices and ensure better
   * distribution of edges across the graph. */
  uint64_t val0, val1;
  {
    mrg_state new_state = state;
    /* Skip ahead 50 * 2^7 iterations to decorrelate from seed */
    mrg_skip(&new_state, 50, 7, 0);
    /* Build val0 as a 64-bit value from two 32-bit RNG outputs */
    val0 = mrg_get_uint_orig(&new_state);
    val0 *= UINT64_C(0xFFFFFFFF);
    val0 += mrg_get_uint_orig(&new_state);
    /* Build val1 similarly for the second scrambling parameter */
    val1 = mrg_get_uint_orig(&new_state);
    val1 *= UINT64_C(0xFFFFFFFF);
    val1 += mrg_get_uint_orig(&new_state);
  }

  /* Parallel loop: generate each edge independently by skipping to its
   * corresponding RNG position. Thread-safe because each iteration uses
   * a distinct RNG state derived from the base state. */
#ifdef _OPENMP
#pragma omp parallel for
#endif
#ifdef __MTA__
#pragma mta assert parallel
#pragma mta block schedule
#endif
  for (ei = start_edge; ei < end_edge; ++ei) {
    /* Create a new RNG state and skip to the position for edge ei */
    mrg_state new_state = state;
    mrg_skip(&new_state, 0, (uint64_t)ei, 0);
    /* Generate a single Kronecker edge and write to output array */
    make_one_edge(nverts, 0, logN, &new_state, edges + (ei - start_edge), val0, val1);
  }
}


/* ====================================
main for basic testing of generation
==================================== */
int main(int argc, char** argv) {
  if (argc != 4) {
    printf("Usage: %s <scale> <edgefactor> <np>\n", argv[0]);
    return 1;
  }

  uint_fast32_t seed[5] = {3, 3, 5, 7, 11};
  int scale = atoi(argv[1]);  /* 2^4 = 16 vertices */
  uint32_t N = pow(2, scale);
  
  int edge_factor = atoi(argv[2]);
  int64_t num_edges = N * edge_factor;
  
  packed_edge* edges = malloc(num_edges*sizeof(packed_edge));
  
  generate_kronecker_range(seed, scale, 0, num_edges, edges);
  
  printf("Generated %ld edges for Kronecker graph with 2^%d vertices:\n", num_edges, scale);
  char filename[256];
  sprintf(filename, "data/raw/kronecker%s.txt", argv[3]); // graph expected to be processed by <np> processes 
  FILE* fp = fopen(filename, "w");
  for (int64_t i = 0; i < num_edges; ++i) {
    fprintf(fp, "%ld %ld\n", edges[i].v0, edges[i].v1);
  }

  free(edges);
  fclose(fp);
  
  return 0;
}