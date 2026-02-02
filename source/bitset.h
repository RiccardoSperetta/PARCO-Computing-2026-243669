#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct {
    uint64_t *bits;
    size_t    n_bits;
    size_t    n_words;
} bitset_t;

bitset_t bitset_create(size_t n_bits);
void bitset_free(bitset_t *bs);
void bitset_clear(bitset_t *bs);
static inline void bitset_set(bitset_t *bs, size_t v) {
    size_t word = v >> 6;          // v / 64 (2^6)
    size_t bit  = v & 63;          // v % 64
    bs->bits[word] |= (1ULL << bit);
}

static inline int bitset_test(const bitset_t *bs, size_t v) {
    size_t word = v >> 6;
    size_t bit  = v & 63;
    return (bs->bits[word] >> bit) & 1ULL;
}