#ifndef HASHSET_H
#define HASHSET_H

#include "cocadautil.hpp"
#include "hashmap.hpp"


/**
 * @file hashset.h
 * @author Paulo Fonseca
 * 
 * @brief Unordered, hashtable-based set - Hashset - ADT.
 */

/**
 * Hashset type
 */
typedef struct _hashset hashset;


/**
 * Hashset iterator type
 */
typedef struct _hashset_iterator hashset_iterator;


/**
 * @brief Constructor.
 * 
 * @param hfunc Hash function pointer.
 * @param eqfunc Equality comparator function pointer.
 */
hashset *hashset_new(hash_func hfunc, equals_func eqfunc);


/**
 * @brief Destructor.
 * @param free_elements Indicates whether stored elements should be disposed.
 */
void hashset_free(hashset *set, bool free_elements);


/**
 * @brief Returns the number of stored elements.
 */
size_t hashset_size();


/**
 * @brief Checks whether the @p set contains a given element @p elt. 
 * @returns true iff @p set contains an element x s.t. 
 *          @p set->eqfunc(x, @p elt) == true
 */
bool hashset_contains(hashset *set, void *elt);


/**
 * @brief Adds an element @p elt to the @p set. 
 */
void hashset_add(hashset *set, void *elt);


/**
 * @brief Removes an element @p elt from the @p set. 
 */
void hashset_remove(hashset *set, void *elt);


/**
 * @brief Returns a new iterator for the @p set.
 */
hashset_iterator *hashset_get_iterator(hashset *set);


/**
 * @brief Iterator destructor.
 * @warning Only the iterator is destroyed. The set is left unmodified.
 */
void hashset_iterator_free(hashset_iterator *it);


/**
 * @brief Indicates whether there are still entries to be iterated over.
 */
bool hashset_iterator_has_next(hashset_iterator *it);


/**
 * @brief Gets the next element of the iteration.
 */
void *hashset_iterator_next(hashset_iterator *it);


#endif
