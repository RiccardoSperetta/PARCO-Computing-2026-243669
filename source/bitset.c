#include "bitset.h"

bitset_t bitset_create(size_t n_bits) {
    bitset_t bs;
    bs.n_bits  = n_bits;
    bs.n_words = (n_bits + 63) / 64; //making sure we don't cut off the last bits 

    bs.bits = calloc(bs.n_words, sizeof(uint64_t));
    if (!bs.bits) {
        perror("calloc bitset");
        exit(EXIT_FAILURE);
    }

    return bs;
}

void bitset_clear(bitset_t *bs) {
    for (int i=0; i<bs->n_words; i++) {
        bs->bits[i] = 0;
    }
}

void bitset_free(bitset_t *bs) {
    free(bs->bits);
    bs->bits = NULL;
}

