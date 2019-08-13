#include "cocadautil.hpp"
#include "hashmap.hpp"
#include "hashset.hpp"

struct _hashset {
    hashmap *hmap;
};

struct _hashset_iterator {
    hashmap_iterator *hmap_it;
};

hashset *hashset_new(size_t (*hash_fn)(void *), bool (*equals_fn)(void *,
                     void *))
{
    hashset *ret;
    ret = NEW(hashset);
    ret->hmap = hashmap_new(hash_fn, equals_fn);
    return ret;
}

void hashset_free(hashset *set, bool free_elements)
{
    if (set->hmap!=NULL) {
        hashmap_free(set->hmap, free_elements, false);
    }
    FREE(set);
}

size_t hashset_size(hashset *set)
{
    return hashmap_size(set->hmap);
}

bool hashset_contains(hashset *set, void *elt)
{
    return hashmap_haskey(set->hmap, elt);
}

void hashset_add(hashset *set, void *elt)
{
    hashmap_set(set->hmap, elt, NULL);
}

void hashset_remove(hashset *set, void *elt)
{
    hashmap_unset(set->hmap, elt);
}


hashset_iterator *hashset_get_iterator(hashset *set)
{
    hashset_iterator *ret;
    ret = NEW(hashset_iterator);
    ret->hmap_it = hashmap_iterator_new(set->hmap);
    return ret;
}

void hashset_iterator_free(hashset_iterator *it)
{
    if (it->hmap_it != NULL) {
        hashmap_iterator_free(it->hmap_it);
    }
    FREE(it);
}

bool hashset_iterator_has_next(hashset_iterator *it)
{
    return hashmap_iterator_has_next(it->hmap_it);
}

void *hashset_iterator_next(hashset_iterator *it)
{
    return hashmap_iterator_next(it->hmap_it)->key;
}

