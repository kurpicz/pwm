#include "bitsandbytes.hpp"
#include "cocadautil.hpp"

void byte_to_str(byte_t b, char *dest)
{
    int i;
    for (i=0; i<BYTESIZE; i++) {
        dest[BYTESIZE-i-1] = '0'+ (b%2) ;
        b >>= 1;
    }
    dest[i]='\0';
}

void byte_reverse(byte_t *b)
{
    *b = (*b<<4) | (*b>>4);
    *b = ((*b<<2) & 0xCC) | ((*b>>2) & 0x33);
    *b = ((*b<<1) & 0xAA) | ((*b>>1) & 0x55);
}


size_t byte_bitcount0(byte_t b)
{
    return byte_bitcount1(~b);
}

size_t byte_bitcount1(byte_t b)
{
    b = ((b>>1) & 0x55)+(b & 0x55);
    b = ((b>>2) & 0x33)+(b & 0x33);
    b = ((b>>4) & 0x0F)+(b & 0x0F);
    return (size_t)b;
}

size_t byte_bitcount(byte_t b, byte_t bit)
{
    if (bit) return byte_bitcount1(b);
    else return byte_bitcount0(b);
}

size_t uint16_bitcount0(uint16_t n)
{
    return uint16_bitcount1(~n);
}

size_t uint16_bitcount1(uint16_t n)
{
    n = ((n>>1) & 0x5555)+(n & 0x5555);
    n = ((n>>2) & 0x3333)+(n & 0x3333);
    n = ((n>>4) & 0x0F0F)+(n & 0x0F0F);
    n = ((n>>8) & 0x00FF)+(n & 0x00FF);
    return (size_t)n;
}

size_t uint16_bitcount(uint16_t n, byte_t bit)
{
    if (bit) 
        return uint16_bitcount1(n);
    else 
        return uint16_bitcount0(n);
}

size_t uint32_bitcount0(uint32_t n)
{
    return uint32_bitcount1(~n);
}

size_t uint32_bitcount1(uint32_t n)
{
    n = ((n>>1) & 0x55555555)+(n & 0x55555555);
    n = ((n>>2) & 0x33333333)+(n & 0x33333333);
    n = ((n>>4) & 0x0F0F0F0F)+(n & 0x0F0F0F0F);
    n = ((n>>8) & 0x00FF00FF)+(n & 0x00FF00FF);
    n = ((n>>16) & 0x0000FFFF)+(n & 0x0000FFFF);
    return (size_t)n;
}

size_t uint32_bitcount(uint32_t n, byte_t bit)
{
    if (bit) 
        return uint32_bitcount1(n);
    else 
        return uint32_bitcount0(n);
}


size_t uint64_bitcount0(uint64_t n)
{
    return uint64_bitcount1(~n);
}


size_t uint64_bitcount1(uint64_t n)
{
    n = ((n>>1) & 0x5555555555555555)+(n & 0x5555555555555555);
    n = ((n>>2) & 0x3333333333333333)+(n & 0x3333333333333333);
    n = ((n>>4) & 0x0F0F0F0F0F0F0F0F)+(n & 0x0F0F0F0F0F0F0F0F);
    n = ((n>>8) & 0x00FF00FF00FF00FF)+(n & 0x00FF00FF00FF00FF);
    n = ((n>>16) & 0x0000FFFF0000FFFF)+(n & 0x0000FFFF0000FFFF);
    n = ((n>>32) & 0x00000000FFFFFFFF)+(n & 0x00000000FFFFFFFF);
    return (size_t)n;
}

size_t uint64_bitcount(uint64_t n, byte_t bit)
{
    if (bit) 
        return uint64_bitcount1(n);
    else 
        return uint64_bitcount0(n);
}


/*
size_t byte_t_select_1(byte_t b, size_t j)
{
    size_t i;
    // revert original number
    byte_t_revert(&b);
    // erase the j-1 last ones
    for (i=1; i<j; i++) {
        b &= (b-1);
    }
    // revert back
    byte_t_revert(&b);
    // count leading zeros
    i = 0;
    if (b==0) return BYTESIZE;
    if (b<=0x0F) {b<<=4; i+=4;}
    if (b<=0x3F) {b<<=2; i+=2;}
    if (b<=0x7F) {i+=1;}
    return i;
}
*/

size_t byte_rank0(byte_t b, size_t pos)
{
    return byte_rank1(~b, pos);
}

size_t byte_rank1(byte_t b, size_t pos)
{
    return byte_bitcount1(b>>(pos<BYTESIZE?BYTESIZE-pos-1:0));
}

size_t byte_rank(byte_t b, size_t pos, byte_t bit)
{
    if (bit) 
        return byte_rank1(b, pos);
    else 
        return byte_rank0(b, pos);
}

size_t byte_select0(byte_t b, size_t rank)
{
    return byte_select1(~b, rank);
}

size_t byte_select1(byte_t b, size_t rank)
{
    size_t i=0;
    if (b==0 || rank==0) 
        return BYTESIZE;
    while (b && rank>=1) {
        // count leading zeros
        if (b <= 0x0F) {
            b<<=4;
            i+=4;
        }
        if (b <= 0x3F) {
            b<<=2;
            i+=2;
        }
        if (b <= 0x7F) {
            b<<=1;
            i+=1;
        }
        rank--;
        b&=0x7F;
    }
    return (rank==0)?i:BYTESIZE;
}

size_t byte_select(byte_t b, size_t rank, byte_t bit)
{
    if (bit) 
        return byte_select1(b, rank);
    else 
        return byte_select0(b, rank);
}
