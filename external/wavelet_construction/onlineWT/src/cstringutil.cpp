#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "cocadautil.hpp"
#include "cstringutil.hpp"

struct _cstrfreader {
    FILE *src;
};


char *cstr_new(size_t len)
{
    char *ret;
    ret = (char *) malloc((len+1)*sizeof(char));
    memset(ret, '\0', (len+1)*sizeof(char));
    return ret;
}

void cstr_fill(char *str, size_t from, size_t to, char c)
{
    for (size_t i=from; i<to; i++) {
        str[i] = c;
    }
}

void cstr_clear(char *str, size_t len)
{
    memset(str, '\0', (len+1)*sizeof(char));
}

char *cstr_trim(char *str, size_t from,  size_t to)
{
    str = (char*)memmove(str, str+from, to-from);
    str = (char*)realloc(str, to-from+1);
    str[to-from] = '\0';
    return str;
}

char *cstr_trim_to_len(char *str, size_t len)
{
    return cstr_trim(str, 0, len);
}


void cstr_revert(char *str, size_t len)
{
    size_t i=0, j=len-1;
    char c;
    while (i<j) {
        c = str[i];
        str[i] = str[j];
        str[j] = c;
        i++;
        j--;
    }
}


static const char *_digits = "0123456789abcdef";

static const size_t _MAX_SIZE_T = ~0x0;


void sizet_to_cstr(char *dest, size_t val, char base)
{
    size_t b, l;
    switch (base) {
    case 'h':
        b = 16;
        l = (size_t)ceil(log2(_MAX_SIZE_T)/4.0);
        break;
    case 'o':
        b = 8;
        l = (size_t)ceil(log2(_MAX_SIZE_T)/3.0);
        break;
    case 'b':
        b = 2;
        l = (size_t)ceil(log2(_MAX_SIZE_T));
        break;
    default:
        b = 10;
        l = (size_t)ceil(log10(_MAX_SIZE_T));
        break;
    }
    memset(dest, '0', l);
    dest[l]='\0';
    while (val) {
        dest[--l] = _digits[val%b];
        val /= b;
    }
}


void *ptr_to_cstr(void *ptr, char *dest)
{
    size_t val, b, l;
    val = (size_t) ptr;
    b = 16;
    l = (size_t)ceil(log2(_MAX_SIZE_T)/4.0);
    if (dest==NULL) {
        dest = cstr_new(l);
    }
    memset(dest, '0', l);
    dest[l]='\0';
    while (val) {
        dest[--l] = _digits[val%b];
        val /= b;
    }
    return dest;
}
