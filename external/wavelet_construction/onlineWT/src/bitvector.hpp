#ifndef BITVECTOR_H
#define BITVECTOR_H

#include <stddef.h>


typedef struct _bitvector bitvector;

bitvector *bitvector_new(size_t capacity);

void bitvector_free(bitvector *bv);

size_t bitvec_len(bitvector *bv);

byte_t bitvec_get_bit (bitvector *bv, size_t pos);

size_t bitvec_count(bitvector *bv, byte_t bit);

void bitvec_append (bitvector *bv, byte_t bit);

void bitvec_append_n (bitvector *bv, size_t n, byte_t bit);

byte_t *bitvec_detach (bitvector *bv);

void bitvec_print_cap(bitvector *bv);

size_t bitvector_cap(bitvector *bv);

void bitvec_print(bitvector *bv, size_t bytes_per_row);


#endif