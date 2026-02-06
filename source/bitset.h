#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct {
    uint64_t *bits;
    size_t    n_bits;
    size_t    n_words;
} bitset_t;

//Creates empty bitset
bitset_t bitset_create(size_t n_bits);

//Deallocates bitset
void bitset_free(bitset_t *bs);

//Sets all bits to 0
void bitset_clear(bitset_t *bs);

//Atomically sets the bit of a bitset to 1
static inline void bitset_set(bitset_t *bs, size_t v) {
    size_t word = v >> 6;          // v / 64 (2^6)
    size_t bit  = v & 63;          // v % 64
    bs->bits[word] |= (1ULL << bit);
}

//Atomically sets the bit of a bitset to 1
static inline void bitset_set_atomic(bitset_t *bs, size_t bit) {
    size_t word = bit / 64;
    size_t offset = bit % 64;
    uint64_t mask = (1ULL << offset);
    
    __atomic_fetch_or(&bs->bits[word], mask, __ATOMIC_RELAXED); //exat same thing of line 25
}

//Returns value of bit indicated by 'v' as if the bitset was an array of bits
static inline int bitset_test(const bitset_t *bs, size_t v) {
    size_t word = v >> 6;
    size_t bit  = v & 63;
    return (bs->bits[word] >> bit) & 1ULL;
}