/*
 * Copyright (C) 2015-  Paulo G.S. da Fonseca
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "arrayutil.hpp"
#include "cocadautil.hpp"
#include "hashmap.hpp"


static size_t DEFAULT_CAPACITY = 17;

static double DEFAULT_MAX_LOAD_FACTOR = 0.5;


/*
struct _hashmap_entry_t
{
    void *key;
    void *val;
};
*/

hashmap_entry *hashmap_entry_new(void *key, void *val)
{
    hashmap_entry *ret;
    ret = NEW(hashmap_entry);
    ret->key = key;
    ret->val = val;
    return ret;
}



struct _hashmap_iterator_t {
    hashmap *source;
    size_t index;
};


struct _hashmap {
    size_t initial_capacity;

    /*
     * The current capacity.
     * When the number of occupied slots equals
     * load_factor\*capacity, the table is expanded and
     * a rehash takes place.
     */
    size_t capacity;

    /*
     * The (logical) size of the table, i.e. the number of
     * stored elements.
     */
    size_t size;

    /*
     * The number of physically occupied slots.
     * A slot may be logically empty but
     * physically occupied  because of a
     * 'soft' deletion which just marks the
     * slot as deleted without necessarily erasing it.
     */
    size_t occupied;

    /*
     * The key comparator function.
     * Any two keys which compare as equal by this function
     * are required to have the same hash value computed
     * by the hash function provided to this table.
     */
    equals_func keyequals;

    /*
     * The hash function which computes an integer 'hash value'
     * for a given key.
     * This function is not strictly required to be one-to-one, although
     * this is usually a desired property.
     * However, any two keys which compare as equal by key comparator
     * function of this hashmap are required to give the
     * same hash value.
     */
    hash_func hash;

    /*
     * The actual array of (key,val) entries.
     */
    hashmap_entry **entries;
};


hashmap *hashmap_new(hash_func hfunc, equals_func eqfunc)
{
    return hashmap_new_with_capacity(hfunc, eqfunc, DEFAULT_CAPACITY);
}


hashmap *hashmap_new_with_capacity(hash_func hfunc, equals_func eqfunc,
                                   size_t capacity)
{
    hashmap *ret;
    size_t i;

    assert(capacity>0);

    ret = NEW(hashmap);

    ret->initial_capacity = capacity;
    ret->capacity = capacity;
    ret->occupied = 0;
    ret->size = 0;
    ret->hash = hfunc;
    ret->keyequals = eqfunc;

    ret->entries = NEW_ARRAY(hashmap_entry *, ret->capacity);
    for (i=0; i<ret->capacity; i++) {
        ret->entries[i] = NULL;
    }
    return ret;
}

void hashmap_free(hashmap *hashmap, bool free_keys, bool free_vals)
{
    size_t i;
    if (free_keys) {
        for (i=0; i<hashmap->capacity; i++) {
            if (hashmap->entries[i]!=NULL) {
                FREE(hashmap->entries[i]->key);
            }
        }
    }
    if (free_vals) {
        for (i=0; i<hashmap->capacity; i++) {
            if (hashmap->entries[i]!=NULL) {
                FREE(hashmap->entries[i]->val);
            }
        }
    }
    for (i=0; i<hashmap->capacity; i++) {
        FREE(hashmap->entries[i]);
    }
    FREE(hashmap->entries);
    FREE(hashmap);
}


size_t _hash(hashmap *hashmap, void *key)
{
    size_t pos;
    pos = hashmap->hash(key) % hashmap->capacity;
    while (hashmap->entries[pos]!=NULL &&
            (hashmap->entries[pos]->key == NULL ||                          // skip soft deleted entries
             !hashmap->keyequals(key,
                                 hashmap->entries[pos]->key))) { // or positions occupied with different keys
        pos = ((pos + 1) % hashmap->capacity);
    }
    return pos;
}

void _rehash(hashmap *hashmap)
{
    size_t newcapacity, k, pos;
    hashmap_entry **newkeys;

    for (newcapacity=hashmap->initial_capacity;
            hashmap->size>=DEFAULT_MAX_LOAD_FACTOR*newcapacity;
            newcapacity=2*newcapacity+1);
    newkeys = NEW_ARRAY(hashmap_entry *, newcapacity);
    for (k=0; k<newcapacity; k++) {
        newkeys[k] = NULL;
    }

    for (k=0; k<hashmap->capacity; k++) {
        if (hashmap->entries[k] == NULL
                || hashmap->entries[k]->key ==
                NULL) // entries with null keys (= deleted) aren't rehashed
            continue;
        pos = hashmap->hash(hashmap->entries[k]->key) % newcapacity;
        while (newkeys[pos]!=NULL
                && !hashmap->keyequals(hashmap->entries[k]->key, newkeys[pos]->key)) {
            pos = ((pos + 1) % newcapacity);
        }
        newkeys[pos] = hashmap->entries[k];
    }

    FREE(hashmap->entries);
    hashmap->capacity = newcapacity;
    hashmap->entries = newkeys;
}


bool hashmap_haskey(hashmap *hashmap, void *key)
{
    return (hashmap->entries[_hash(hashmap, key)] != NULL);
}


void *hashmap_get(hashmap *hashmap, void *key)
{
    size_t pos;
    assert(key != NULL);
    pos = _hash(hashmap, key);
    return ((hashmap->entries[pos]!=NULL)?hashmap->entries[pos]->val:NULL);
}

void hashmap_set(hashmap *hashmap, void *key, void *val)
{
    assert(key != NULL);
    size_t pos;
    if (hashmap->occupied >= (DEFAULT_MAX_LOAD_FACTOR * hashmap->capacity)) {
        _rehash(hashmap);
    }
    pos = _hash(hashmap, key);

    if (hashmap->entries[pos]==NULL) {
        hashmap->entries[pos] = NEW(hashmap_entry);
        hashmap->occupied++;
        hashmap->size++;
    }
    hashmap->entries[pos]->key = key;
    hashmap->entries[pos]->val = val;
}


void hashmap_unset(hashmap *hashmap, void *key)
{
    size_t pos;
    assert(key != NULL);
    pos = _hash(hashmap, key);
    if (hashmap->entries[pos]!=NULL) {
        hashmap->entries[pos]->key = NULL;
        hashmap->entries[pos]->val = NULL;
        hashmap->size--;
    }
}

size_t hashmap_size(hashmap *map)
{
    return map->size;
}



hashmap_iterator *hashmap_iterator_new(hashmap *source)
{
    hashmap_iterator *ret;
    ret = NEW(hashmap_iterator);
    ret->source = source;
    ret->index = 0;
    while ((ret->index < ret->source->capacity) &&
            (ret->source->entries[ret->index] == NULL ||
             ret->source->entries[ret->index]->key == NULL)) {
        ret->index++;
    }
    return ret;
}


void hashmap_iterator_free(hashmap_iterator *it)
{
    FREE(it);
}

bool hashmap_iterator_has_next(hashmap_iterator *it)
{
    return it->index < it->source->capacity;
}

hashmap_entry *hashmap_iterator_next(hashmap_iterator *it)
{
    hashmap_entry *ret;
    if (it->index >= it->source->capacity) {
        return NULL;
    }
    ret = hashmap_entry_new(it->source->entries[it->index]->key,
                            it->source->entries[it->index]->val);
    do {
        it->index++;
    }
    while ((it->index < it->source->capacity) &&
            (it->source->entries[it->index] == NULL ||
             it->source->entries[it->index]->key == NULL));
    return ret;
}

