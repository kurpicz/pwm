#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "arrayutil.hpp"
#include "bitarray.hpp"
#include "bitsandbytes.hpp"
#include "bitvector.hpp"
#include "cocadautil.hpp"
#include "mathutil.hpp"

#define GROW_BY  1.5f
#define BYTE_MSB  0x80


struct _bitvector {
    byte_t *bits;
    size_t len;
    size_t cap;
    size_t byte_cap;
    //size_t count1;
};

bitvector *bitvector_new(size_t capacity)
{
    bitvector *bv = NEW(bitvector);
    bv->len = 0;
    bv->byte_cap = MAX(16, (size_t)multceil(capacity, BYTESIZE));
    bv->cap = bv->byte_cap*BYTESIZE;
    bv->bits = (byte_t*)malloc(bv->byte_cap);
    memset(bv->bits, 0x00, bv->byte_cap);
    //bv->count1 = 0;
    return bv;
}


void bitvector_free(bitvector *bv) 
{
    FREE(bv->bits);
    FREE(bv);
}


size_t bitvec_len(bitvector *bv) {
    return bv->len;
}


byte_t bitvec_get_bit (bitvector *bv, size_t pos)
{
    byte_t bit_pos = pos%BYTESIZE;
    return ((bv->bits[pos/BYTESIZE] & (BYTE_MSB>>bit_pos)) >> (BYTESIZE-bit_pos-1));
}


static inline size_t _bitvec_count1(bitvector *bv)
{
    size_t byte_pos = 0;
    size_t used_bytes = (size_t)multceil(bv->len, BYTESIZE);
    size_t cnt = 0;
#if BYTEWORDSIZE==8
    while (byte_pos+8 < used_bytes) {
        cnt += uint64_bitcount1(*((uint64_t *)(bv->bits+byte_pos)));
        byte_pos += 8;
    }
    while (byte_pos+4 < used_bytes) {
        cnt += uint32_bitcount1(*((uint32_t *)(bv->bits+byte_pos)));
        byte_pos += 4;
    }
    while (byte_pos+2 < used_bytes)
    {
        cnt += uint16_bitcount1(*((uint16_t *)(bv->bits+byte_pos)));
        byte_pos += 2;
    }
#elif BYTEWORDSIZE==4
    while (byte_pos+4 < used_bytes) {
        cnt += uint32_bitcount1(*((uint32_t *)(bv->bits+byte_pos)));
        byte_pos += 4;
    }
    while (byte_pos+2 < used_bytes) {
        cnt += uint16_bitcount1(*((uint16_t *)(bv->bits+byte_pos)));
        byte_pos += 2;
    }
#endif
    while (byte_pos < used_bytes) {
        cnt += byte_bitcount1(bv->bits[byte_pos]);
        byte_pos++;
    }
    return cnt;
}


static inline size_t _bitvec_count0(bitvector *bv) {
    return bv->len - _bitvec_count1(bv);
}

/*
static inline size_t __bitvec_count1(bitvector *bv) {
    return bv->count1;
}

static inline size_t __bitvec_count0(bitvector *bv) {
    return bv->len - bv->count1;
}
*/

typedef size_t (*_bv_cnt_func)(bitvector *);


static _bv_cnt_func _bitvec_count_func[2] = {_bitvec_count0, _bitvec_count1};


size_t bitvec_count(bitvector *bv, byte_t bit)
{
    return _bitvec_count_func[bit](bv);
}


static void _growto(bitvector *bv, size_t new_byte_cap) {
    bv->bits = (byte_t*)realloc(bv->bits, new_byte_cap);
    memset(bv->bits+bv->byte_cap, 0x00, (new_byte_cap - bv->byte_cap));
    bv->byte_cap = new_byte_cap;
    bv->cap = bv->byte_cap*BYTESIZE;
}


void bitvec_append (bitvector *bv, byte_t bit)
{
    bv->bits[bv->len/BYTESIZE] ^= ( ((-bit)^(bv->bits[bv->len/BYTESIZE]))
                                     & (BYTE_MSB>>(bv->len % BYTESIZE)) );
    if ((++bv->len)==bv->cap) 
        _growto(bv, GROW_BY*bv->byte_cap);
}


void bitvec_append_n (bitvector *bv, size_t nbits, byte_t bit)
{
    if (bv->len+nbits >= bv->cap) { // HAS to be >=, not > 
        _growto(bv, GROW_BY*(size_t)(multceil(bv->len+nbits, BYTESIZE)));
    }
    if (bit) {
        size_t nleft=nbits;
        byte_t nxt_bit = bv->len % BYTESIZE;
        byte_t *last_byte = bv->bits + (bv->len/BYTESIZE); 
        size_t m = MIN(nleft, BYTESIZE-nxt_bit);
        if (m==BYTESIZE)
            *(last_byte) = BYTE_MAX;
        else 
            *(last_byte) |= (~(BYTE_MAX<<m)) << (BYTESIZE - nxt_bit - m);  
        last_byte += ((nxt_bit+m)/BYTESIZE);
        nxt_bit = ((nxt_bit+m)%BYTESIZE);
        nleft -= m;
        m = nleft/BYTESIZE;
        memset(last_byte, 0xFF, m);
        last_byte += m;
        nleft -= (m*BYTESIZE);
        *(last_byte) |= ~(BYTE_MAX>>nleft);  
        bv->len += nbits;
        //bv->count1 += nbits;
    }
    else    
        bv->len += nbits;

}


byte_t *bitvec_detach (bitvector *bv)
{
    byte_t *ret = (byte_t *)bv->bits;
    ret = (byte_t*)realloc(ret, (size_t)multceil(bv->len, BYTESIZE));
    FREE(bv);
    return ret;
}

void bitvec_print_cap(bitvector *bv) {
    printf("Capacity: %zu\n", bv->cap);
}

size_t bitvector_cap(bitvector *bv) {
    return bv->cap;
}

void bitvec_print(bitvector *bv, size_t bytes_per_row)
{
    printf("bitvector@%p\n", bv);
    printf("  len     : %zu\n", bv->len);
    printf("  capacity: %zu\n", bv->cap);
    //printf("  count1: %zu\n", bv->count1);
    printf("  data:\n");
    bitarr_print(bv->bits, bv->cap, bytes_per_row);
    printf("end of bitvector@%p\n", bv);
    
}

