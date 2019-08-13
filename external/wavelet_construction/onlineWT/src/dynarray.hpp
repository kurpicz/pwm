#ifndef DYNARRAY_H
#define DYNARRAY_H

#include "bitsandbytes.hpp"
#include "cocadautil.hpp"

typedef struct _dynarray dynarray;

/**
 * @file dynarray.h
 * @author Paulo Fonseca
 *
 * @brief Dynamic array ADT (a.k.a. Vector)
 *        The array can be used to store pointers to external objects or 
 *        values of primitive types such as int, via specific read/write 
 *        functions. 
 */


/**
 * @brief Array constructor.
 * @param typesize The size of the elements to be stored (in bytes).
 */
dynarray *dynarray_new(size_t typesize);


/**
 * @brief Array constructor.
 * @param typesize The size of the elements to be stored (in bytes).
 * @param init_capacity The initial capacity (in # of elements).
 */
dynarray *dynarray_new_with_capacity(size_t typesize, size_t init_capacity);


/**
 * @brief Destructor.
 * @param free_elements Indicates whether referenced objects should be freed.
 */
void dyn_array_free(dynarray *da, bool free_elements);


/**
 * @brief Returns the length, i.e. the # of elements "logically" stored.
 */
size_t darr_len(dynarray *da);


/**
 * @brief Returns the individual size of stored elements (in bytes). 
 */
size_t darr_typesize(dynarray *da);


/**
 * @brief Removes all elements <b>without</b> destroying them.
 */
void darr_clear(dynarray *da);


/**
 * @brief Returns a read-only reference to the current internal array.
 *
 * @warning The internal array can change between calls and the returned 
 *          reference can become NULL or invalid.
 */
void *darr_as_array(dynarray *da);


/**
 * @brief Detaches and returns the current internal array 
*         (including trailing unused positions).
 * 
 * @see darr_trim
 * @warning After this operation, the dynamic array @p da is destroyed.
 */
void *darr_detach(dynarray *da);


/**
 * @brief Returns the element at position @p pos.
 */
void *darr_get(dynarray *da, size_t pos);


/**
 * @brief Sets (overwrites) the element at position @p pos to @p val.
 */
void  darr_set(dynarray *da, size_t pos, void *val);


/**
 * @brief Appends a new element @p val.
 */
void  darr_app(dynarray *da, void *val);


/**
 * @brief Inserts a new element @p val at position @p pos.
 */
void  darr_ins(dynarray *da, size_t pos, void *val);


/**
 * @brief Removes and returns new the element at position @p pos.
 */
void *darr_del(dynarray *da, size_t pos);


int   darr_get_int(dynarray *da, size_t pos);
void  darr_set_int(dynarray *da, size_t pos, int val);
void  darr_app_int(dynarray *da, int val);
void  darr_ins_int(dynarray *da, size_t pos, int val);
int   darr_del_int(dynarray *da, size_t pos);


size_t darr_get_sizet(dynarray *da, size_t pos);
void   darr_set_sizet(dynarray *da, size_t pos, size_t val);
void   darr_app_sizet(dynarray *da, size_t val);
void   darr_ins_sizet(dynarray *da, size_t pos, size_t val);
size_t darr_del_sizet(dynarray *da, size_t pos);


byte_t darr_get_byte(dynarray *da, size_t pos);
void   darr_set_byte(dynarray *da, size_t pos, byte_t val);
void   darr_app_byte(dynarray *da, byte_t val);
void   darr_ins_byte(dynarray *da, size_t pos, byte_t val);
byte_t darr_del_byte(dynarray *da, size_t pos);


/**
 * @brief Radix sort on a generic dynamic array.
 *
 * We assume that each element of the array can be associated to a
 * numeric key vector of key_size positions K = (K[key_size-1],...,K[0]),
 * where each key position K[j] assumes one of max_key integer values in
 * the range 0 <= K[j]< max_key.
 * Then this function performs a radix sort on the elements of the array based
 * on their key vectors, with K[0] being the least significant position, and
 * K[key_size-1] the most significant. That is the elements end up sorted by
 * K[key_size-1], then K[key_size-2], and so forth, downto K[0].
 *
 * @param da The array to sort
 * @param key_fn A pointer to a function that computes K[j] from a given element
 * @param key_size The size of the key vector
 * @param max_key The noninclusive maximum value for each key position
 */
void darr_radixsort(dynarray *da, size_t (*key_fn)(void *, size_t),
                    size_t key_size, size_t max_key);
                    
                    
#endif
