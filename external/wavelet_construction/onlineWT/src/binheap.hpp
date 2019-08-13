#ifndef BINHEAP_H
#define BINHEAP_H

#include <stdlib.h>

#include "cocadautil.hpp"

/**
 * @file binheap.h
 * @author Paulo Fonseca
 *
 * @brief Binary heap.
 *
 * This generic binary heap maintains a dynamic collection of items
 * over an ordered set. It supports at least insertion and min/max
 * extraction in O(log(n)) time.
 *
 * This heap will usually store only references to actual "objects".
 * However, it can also store copies of primitive types like numbers
 * and characters through type-specific push/pop operations.
 */

/**
 * @brief The mode of the heap: in a MIN_HEAP (MAX_HEAP) the elements are
 * extracted in MIN-first (resp. MAX-first) order.
 */
typedef enum {
    MIN_HEAP = 0,
    MAX_HEAP = 1
} heap_mode;


/**
 * @brief Binary heap type.
 */
typedef struct _binheap binheap;


/**
 * @brief Creates a new empty binary heap.
 * @param typesize The size (in bytes) of the individual elements.
 * @param mode The mode of the heap
 */
binheap *binheap_new(int (*comp_fn)(void *, void *), size_t typesize,
                     heap_mode mode);


/**
 * @brief Destructor.
 * @param free_elements Indicates whether stored elements should be
 * individually freed.
 */
void binheap_free(binheap *heap, bool free_elements);


/**
 * @brief Returns the number of elements stored in the heap.
 */
size_t binheap_size(binheap *heap);


/**
 * @brief Stores a new element in the heap.
 * @param elt A pointer to the element to be stored. Only this reference is
 * actually stored. Upon heap destruction, the pointed memory is only freed
 * if free_elements is set to 1.
 *
 * @see binheap_free
 */
void binheap_push(binheap *heap, void *elt);


/**
 * @brief Removes and returns the reference MIN/MAX element stored in
 *        the heap depending on the heap mode.
 *
 * @see heap_mode
 */
void *binheap_pop(binheap *heap);


/**
 * @brief Pushes an int onto the heap.
 */
void binheap_push_int(binheap *heap, int val);


/**
 * @brief Pops the MIN/MAX int from the heap.
 */
int binheap_pop_int(binheap *heap);


/**
 * @brief Pushes a size_t onto the heap.
 */
void binheap_push_size(binheap *heap, size_t val);


/**
 * @brief Pops the MIN/MAX size_t from the heap.
 */
size_t binheap_pop_size(binheap *heap);


#endif
