#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "arrayutil.hpp"
#include "cocadautil.hpp"
#include "dynarray.hpp"


const static size_t DYNARRAY_DEFAULT_CAPACITY = 4;

typedef struct _dynarray {
    size_t typesize;
    size_t len;
    size_t capacity;
    void *data;
}
dynarray;

dynarray *dynarray_new(size_t typesize)
{
    return dynarray_new_with_capacity(typesize, DYNARRAY_DEFAULT_CAPACITY);
}

dynarray *dynarray_new_with_capacity(size_t typesize, size_t init_capacity)
{
    dynarray *ret;
    ret = NEW(dynarray);
    ret->typesize = typesize;
    ret->capacity = init_capacity;
    ret->len = 0;
    ret->data = NEW_ARRAY(void, ret->capacity*ret->typesize);
    return ret;

}

static void _double(dynarray *da)
{
    da->capacity = MAX(2*da->capacity, 1);
    da->data = realloc(da->data, da->capacity*da->typesize);
}

void dyn_array_free(dynarray *da, bool free_elements)
{
    if (free_elements) {
        for (size_t i=0; i<da->len; i++) {
            FREE(((void **)da->data)[i]);
        }
    }
    FREE(da->data);
    FREE(da);
}

size_t darr_len(dynarray *da)
{
    return da->len;
}

size_t darr_typesize(dynarray *da)
{
    return da->typesize;
}

void darr_clear(dynarray *da)
{
    da->len = 0;
}

void *darr_as_array(dynarray *da)
{
    return da->data;
}


void *darr_detach(dynarray *da)
{
    void *data = da->data;
    FREE(da);
    return data;
}


void *darr_get(dynarray *da, size_t pos)
{
    return ((void **)da->data)[pos];
}

void darr_set(dynarray *da, size_t pos, void *val)
{
    ((void **)da->data)[pos] = val;
}

void darr_app(dynarray *da, void *val)
{
    if (da->len==da->capacity) {
        _double(da);
    }
    ((void **)da->data)[da->len++] = val;
}

void darr_ins(dynarray *da, size_t pos, void *val)
{
    if (da->len==da->capacity) {
        _double(da);
    }
    for (size_t i=da->len; i>pos; i--) {
        ((void **)da->data)[i] = ((void **)da->data)[i-1];
    }
    ((void **)da->data)[pos] = val;
    da->len++;
}

void *darr_del(dynarray *da, size_t pos)
{
    void *ret = ((void **)da->data)[pos];
    for (size_t i=pos; i<da->len-1; i++) {
        ((void **)da->data)[i] = ((void **)da->data)[i+1];
    }
    da->len--;
    return ret;
}

int darr_get_int(dynarray *da, size_t pos)
{
    return ((int *)da->data)[pos];
}

void darr_set_int(dynarray *da, size_t pos, int val)
{
    ((int *)da->data)[pos] = val;
}

void darr_app_int(dynarray *da, int val)
{
    if (da->len==da->capacity) {
        _double(da);
    }
    ((int *)da->data)[da->len++] = val;
}

void darr_ins_int(dynarray *da, size_t pos, int val)
{
    if (da->len==da->capacity) {
        _double(da);
    }
    for (size_t i=da->len; i>pos; i--) {
        ((int *)da->data)[i] = ((int *)da->data)[i-1];
    }
    ((int *)da->data)[pos] = val;
    da->len++;
}

int darr_del_int(dynarray *da, size_t pos)
{
    int ret = ((int *)da->data)[pos];
    if (da->len==da->capacity) {
        _double(da);
    }
    for (size_t i=pos; i<da->len-1; i++) {
        ((int *)da->data)[i] = ((int *)da->data)[i+1];
    }
    da->len--;
    return ret;
}


size_t darr_get_sizet(dynarray *da, size_t pos)
{
    return ((size_t *)da->data)[pos];
}

void darr_set_sizet(dynarray *da, size_t pos, size_t val)
{
    ((size_t *)da->data)[pos] = val;
}

void darr_app_sizet(dynarray *da, size_t val)
{
    if (da->len==da->capacity) {
        _double(da);
    }
    ((size_t *)da->data)[da->len++] = val;
}

void darr_ins_sizet(dynarray *da, size_t pos, size_t val)
{
    if (da->len==da->capacity) {
        _double(da);
    }
    for (size_t i=da->len; i>pos; i--) {
        ((size_t *)da->data)[i] = ((size_t *)da->data)[i-1];
    }
    ((size_t *)da->data)[pos] = val;
    da->len++;
}

size_t darr_del_sizet(dynarray *da, size_t pos)
{
    size_t ret = ((size_t *)da->data)[pos];
    for (size_t i=pos; i<da->len-1; i++) {
        ((size_t *)da->data)[i] = ((size_t *)da->data)[i+1];
    }
    da->len--;
    return ret;
}


byte_t darr_get_byte(dynarray *da, size_t pos)
{
    return ((byte_t *)da->data)[pos];
}

void darr_set_byte(dynarray *da, size_t pos, byte_t val)
{
    ((byte_t *)da->data)[pos] = val;
}

void darr_app_byte(dynarray *da, byte_t val)
{
    if (da->len==da->capacity) {
        _double(da);
    }
    ((byte_t *)da->data)[da->len++] = val;
}

void darr_ins_byte(dynarray *da, size_t pos, byte_t val)
{
    if (da->len==da->capacity) {
        _double(da);
    }
    for (size_t i=da->len; i>pos; i--) {
        ((byte_t *)da->data)[i] = ((byte_t *)da->data)[i-1];
    }
    ((byte_t *)da->data)[pos] = val;
    da->len++;
}

byte_t darr_del_byte(dynarray *da, size_t pos)
{
    byte_t ret = ((byte_t *)da->data)[pos];
    if (da->len==da->capacity) {
        _double(da);
    }
    for (size_t i=pos; i<da->len-1; i++) {
        ((byte_t *)da->data)[i] = ((byte_t *)da->data)[i+1];
    }
    da->len--;
    return ret;
}





void darr_radixsort(dynarray *da, size_t (*key_fn)(void *, size_t),
                    size_t key_size, size_t max_key)
{
    size_t i, d, k, n;
    void **dacpy;
    size_t *count;
    n = darr_len(da);
    dacpy = NEW_ARRAY(void *, n);
    count = NEW_ARRAY(size_t, max_key);
    for (d=0; d<key_size; d++) {
        for (k=0; k<max_key; k++) {
            count[k] = 0;
        }
        for (i=0; i<n; i++) {
            k = key_fn(darr_get(da, i), d);
            count[k]++;
        }
        for (k=1; k<max_key; k++) {
            count[k] += count[k-1];
        }
        for (i=n; i>0; i--) {
            k = key_fn(darr_get(da, i-1), d);
            dacpy[--count[k]] = darr_get(da, i-1);
        }
        for (i=0; i<n; i++) {
            darr_set(da, i, dacpy[i]);
        }
    }
    FREE(dacpy);
    FREE(count);
}
