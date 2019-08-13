#include <string.h>

#include "arrayutil.hpp"
#include "binheap.hpp"
#include "cocadautil.hpp"

static const size_t DEFAULT_HEAP_CAPACITY = 4;

struct _binheap {
    size_t typesize;
    size_t size;
    size_t capacity;
    heap_mode mode;
    int (*comp_fn)(void *, void *);
    void *data;
};

binheap *binheap_new(int (*comp_fn)(void *, void *), size_t typesize,
                     heap_mode mode)
{
    binheap *heap;
    heap = NEW(binheap);
    heap->typesize = typesize;
    heap->size = 0;
    heap->capacity = DEFAULT_HEAP_CAPACITY;
    heap->mode = mode;
    heap->comp_fn = comp_fn;
    heap->data = NEW_ARRAY(void, heap->capacity*heap->typesize);
    return heap;
}

void binheap_free(binheap *heap, bool free_elements)
{
    if (free_elements) {
        for (size_t i=0; i<heap->size; i++) {
            FREE(((void **)heap->data)[i]);
        }
    }
    FREE(heap->data);
    FREE(heap);
}


size_t binheap_size(binheap *heap)
{
    return heap->size;
}

static void _double(binheap *heap)
{
    heap->capacity = MAX(2*heap->capacity, 1);
    heap->data = realloc(heap->data, heap->capacity*heap->typesize);
}


static size_t _bubble_up(binheap *heap, size_t pos)
{
    size_t i=pos;
    void *par, *cur, *swp;
    swp = malloc(heap->typesize);
    switch (heap->mode) {
    case MIN_HEAP:
        while ( i>0 &&
                heap->comp_fn(cur=heap->data+(i*heap->typesize) ,
                              par=heap->data+(((i-1)/2)*heap->typesize))<0 ) {
            memcpy(swp, cur, heap->typesize);
            memcpy(cur, par, heap->typesize);
            memcpy(par, swp, heap->typesize);
            i = (i-1)/2;
        }
        break;
    case MAX_HEAP:
        while ( i>0 &&
                heap->comp_fn(cur=heap->data+(i*heap->typesize) ,
                              par=heap->data+(((i-1)/2)*heap->typesize))>0 ) {
            memcpy(swp, cur, heap->typesize);
            memcpy(cur, par, heap->typesize);
            memcpy(par, swp, heap->typesize);
            i = (i-1)/2;
        }
        break;
    }
    free(swp);
    return i;
}


static size_t _bubble_down(binheap *heap, size_t pos)
{
    size_t i, l, r, m, n;
    void *swp;
    swp = malloc(heap->typesize);
    n = binheap_size(heap);
    i = pos;
    while (true) {
        m = i;
        l = (2*i)+1;
        r = (2*i)+2;
        switch (heap->mode) {
        case MIN_HEAP:
            if ( l < n &&
                    heap->comp_fn( heap->data+(l*heap->typesize),
                                   heap->data+(m*heap->typesize) ) < 0 ) {
                m = l;
            }
            if ( r < n &&
                    heap->comp_fn( heap->data+(r*heap->typesize),
                                   heap->data+(m*heap->typesize) ) < 0 ) {
                m = r;
            }
            break;
        case MAX_HEAP:
            if ( l < n &&
                    heap->comp_fn( heap->data+(l*heap->typesize),
                                   heap->data+(m*heap->typesize) ) > 0 ) {
                m = l;
            }
            if ( r < n &&
                    heap->comp_fn( heap->data+(r*heap->typesize),
                                   heap->data+(m*heap->typesize) ) > 0 ) {
                m = r;
            }
            break;
        }
        if ( m != i ) {
            memcpy(swp, heap->data+(i*heap->typesize), heap->typesize);
            memcpy(heap->data+(i*heap->typesize), heap->data+(m*heap->typesize),
                   heap->typesize);
            memcpy(heap->data+(m*heap->typesize), swp, heap->typesize);
            i = m;
        }
        else {
            break;
        }
    }
    free(swp);
    return i;
}


void binheap_push(binheap *heap, void *elt)
{
    if ( heap->size == heap->capacity ) {
        _double(heap);
    }
    ((void **)heap->data)[heap->size] = elt;
    heap->size++;
    _bubble_up(heap, heap->size-1);
}


void *binheap_pop(binheap *heap)
{
    void *first;
    first = ((void **)heap->data)[0];
    heap->size--;
    if (heap->size > 0) {
        ((void **)heap->data)[0]    = ((void **)heap->data)[heap->size];
        _bubble_down(heap, 0);
    }
    return first;
}

void binheap_push_int(binheap *heap, int val)
{
    if ( heap->size == heap->capacity ) {
        _double(heap);
    }
    ((int *)heap->data)[heap->size] = val;
    heap->size++;
    _bubble_up(heap, heap->size-1);
}

int binheap_pop_int(binheap *heap)
{
    int first;
    first = ((int *)heap->data)[0];
    heap->size--;
    if (heap->size > 0) {
        ((int *)heap->data)[0]    = ((int *)heap->data)[heap->size];
        _bubble_down(heap, 0);
    }
    return first;
}


void binheap_push_size(binheap *heap, size_t val)
{
    if ( heap->size == heap->capacity ) {
        _double(heap);
    }
    ((size_t *)heap->data)[heap->size] = val;
    heap->size++;
    _bubble_up(heap, heap->size-1);
}


size_t binheap_pop_size(binheap *heap)
{
    size_t first;
    first = ((size_t *)heap->data)[0];
    heap->size--;
    if (heap->size > 0) {
        ((size_t *)heap->data)[0] = ((size_t *)heap->data)[heap->size];
        _bubble_down(heap, 0);
    }
    return first;
}
