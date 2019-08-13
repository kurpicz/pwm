#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "alphabet.hpp"
#include "arrayutil.hpp"
#include "cocadautil.hpp"
#include "cstringutil.hpp"
#include "hashmap.hpp"
//#include "marshall.h"
#include "strstream.hpp"


typedef enum {
    _AB_PLAIN = 0,
    _AB_MAP   = 1,
    _AB_FUNC  = 2
} _ab_type;


struct _alphabet {
    byte_t mode;
    size_t size;
    char   *letters;
    union {
        size_t *arr;
        char_rank_func func;
        hashmap *map;
    } ranks;
};


alphabet *alphabet_new(const size_t size, const char *letters)
{
    alphabet *ret;
    ret =  NEW(alphabet);
    ret->mode = _AB_PLAIN;
    ret->size = size;
    ret->letters = cstr_new(size);
    strncpy(ret->letters, letters, size);
    ret->ranks.arr = NEW_ARRAY(size_t, UCHAR_MAX);
    FILL_ARRAY(ret->ranks.arr, 0, UCHAR_MAX, size);
    for (size_t i=0; i<size; i++)
        ret->ranks.arr[(size_t)((unsigned char)letters[i])]=i;
    return ret;
}


alphabet *alphabet_new_with_rank_func( const size_t size, const char *letters,
                                       char_rank_func crfunc )
{
    alphabet *ret;
    ret =  NEW(alphabet);
    ret->mode = _AB_FUNC;
    ret->size = size;
    ret->letters = cstr_new(size);
    strncpy(ret->letters, letters, size);
    ret->ranks.func = crfunc;
    return ret;
}


static size_t _char_hash(void *c)
{
    return (size_t)(*((char *)c));
}


bool _char_equals(void *c1, void *c2)
{
    return ((*((char *)c1))==(*((char *)c2)));
}


alphabet *alphabet_new_with_rank_map(const size_t size, const char *letters)
{
    alphabet *ret;
    ret =  NEW(alphabet);
    ret->mode = _AB_MAP;
    ret->size = size;
    ret->letters = cstr_new(size);
    strncpy(ret->letters, letters, size);
    ret->ranks.map = hashmap_new_with_capacity( &_char_hash, &_char_equals,
                                                2*size );
    for (size_t i=0; i<ret->size; i++) {
        size_t *rk = NEW(size_t);
        *rk = i;
        hashmap_set(ret->ranks.map, ret->letters+i, rk);
    }
    return ret;
}


void ab_print(alphabet *ab)
{
    printf("alphabet@%p\n",(void *)ab);
    printf("  size = %zu\n", ab->size);
    printf("  letters = %s\n", ab->letters);
}


void alphabet_free(alphabet *ab)
{
    switch (ab->mode) {
    case _AB_PLAIN:
        FREE(ab->ranks.arr);
        break;
    case _AB_MAP:
        hashmap_free(ab->ranks.map, false, true);
        break;
    }
    FREE(ab->letters);
    FREE(ab);
}


size_t ab_size(alphabet *ab)
{
    return ab->size;
}


bool ab_contains(alphabet *ab, char c)
{
    return (bool)(ab_rank(ab, c) < ab->size);
}



char ab_char(alphabet *ab, size_t index)
{
    return ab->letters[index];
}


size_t ab_rank(alphabet *ab, char c)
{
    switch (ab->mode) {
    case _AB_PLAIN:
        return ab->ranks.arr[(unsigned char)c];
        break;
    case _AB_FUNC:
        return ab->ranks.func(c);
        break;
    case _AB_MAP:
        if (! hashmap_haskey(ab->ranks.map, &c) )
            return ab->size;
        return *((size_t *)(hashmap_get(ab->ranks.map, &c)));
        break;
    }
    return ab->size;
}


size_t *ab_count(alphabet *ab, char *str, size_t slen)
{
    size_t absz = ab_size(ab);
    size_t *counts = NEW_ARRAY(size_t, absz);
    FILL_ARRAY(counts, 0, absz, 0);
    for (size_t i=0; i<slen; counts[ab_rank(ab, str[i++])]++);
    return counts;
}


size_t  *ab_cumul_count(alphabet *ab, char *str, size_t slen)
{
    size_t absz = ab_size(ab);
    size_t *cumul_freqs;
    cumul_freqs = NEW_ARRAY(size_t, absz+1);
    FILL_ARRAY(cumul_freqs, 0, absz+1, 0);
    for (size_t i=0; i<slen; i++) 
        cumul_freqs[ab_rank(ab, str[i])+1]++;
    for (size_t i=1; i<=absz; i++)
        cumul_freqs[i] += cumul_freqs[i-1];
    return cumul_freqs;
}


size_t *ab_count_stream(alphabet *ab, strstream *sst)
{
    size_t absz = ab_size(ab);
    size_t *counts = NEW_ARRAY(size_t, absz);
    FILL_ARRAY(counts, 0, absz, 0);
    strstream_reset(sst);
    for (int c; (c=strstream_getc(sst))!=EOF; counts[ab_rank(ab, (char)c)]++);
    return counts;
}


size_t  *ab_cumul_count_stream(alphabet *ab, strstream *sst)
{
    size_t absz = ab_size(ab);
    size_t *cumul_freqs;
    cumul_freqs = NEW_ARRAY(size_t, absz+1);
    FILL_ARRAY(cumul_freqs, 0, absz+1, 0);
    strstream_reset(sst);
    for ( int c; (c=strstream_getc(sst))!=EOF;
            cumul_freqs[ab_rank(ab, (char)c)+1]++ );
    for (size_t i=1; i<=absz; i++)
        cumul_freqs[i] += cumul_freqs[i-1];
    return cumul_freqs;
}


const char *ab_as_str(alphabet *ab)
{
    return ab->letters;
}

/*
void ab_marshall(marshallctx *ctx, alphabet *ab)
{
    marshall_sizet(ctx, (size_t)ab);
    if (marshallctx_has(ctx, ab))
        return;
    marshall_byte(ctx, ab->mode);
    marshall_sizet(ctx, ab->size);
    marshall_str(ctx, ab->letters, ab->size+1);
    if (ab->mode==_AB_FUNC)
        marshall_func(ctx, ab->ranks.func);
    marshallctx_set(ctx, ab, ab);
}


alphabet *ab_unmarshall(marshallctx *ctx)
{
    size_t *uid = NEW(size_t);
    *uid = unmarshall_sizet(ctx);
    if (marshallctx_has(ctx, uid))
        return marshallctx_get(ctx, uid);
    byte_t mode = unmarshall_byte(ctx);
    size_t size = unmarshall_sizet(ctx);
    char *letters = unmarshall_str(ctx);
    char_rank_func crf;
    alphabet *ab = NULL;
    switch (mode) {
    case _AB_PLAIN:
        ab = alphabet_new(size, letters);
        break;
    case _AB_MAP:
        ab = alphabet_new_with_rank_map(size, letters);
        break;
    case _AB_FUNC:
        crf = unmarshall_func(ctx);
        ab = alphabet_new_with_rank_func(size, letters, crf);
        break;
    }
    marshallctx_set(ctx, uid, ab);
    return ab;
}
*/
