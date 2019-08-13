#ifndef ALPHABET_H
#define ALPHABET_H

#include <stddef.h>
#include <stdlib.h>

//#include "marshall.h"
#include "strstream.hpp"


/**
 * @file alphabet.h
 * @author Paulo Fonseca
 *
 * @brief Finite ordered alphabet ADT.
 *
 * An alphabet is a finite ordered set of letters
 * A=(a[0], a[1], ..., a[l-1]) s.t. rank(a[j])==j, and so
 * a[0] < a[1] < ... < a[l-1] where < stands for the lexicographic order.
 */

/**
 * @brief alphabet type
 */
typedef struct _alphabet alphabet;


/**
 * @brief char rank function type
 */
typedef size_t (*char_rank_func)(char c);


/**
 * @brief Creates an alphabet from a string with letters in lexicographic order.
 * Letter rank computation is unspecified.
 * @param size Number of letters.
 * @param letters String with letters in lexicographic order.
 */
alphabet *alphabet_new(const size_t size, const char *letters);


/**
 * @brief Creates an alphabet from a string with letters in lexicographic order.
 *        Letter rank computation will be done via the provided function.
 * @param size Number of letters.
 * @param letters String with letters in lexicographic order.
 * @param rfunc Pointer to a rank function that computes the position of
 * a letter in the alphabet.
 */
alphabet *alphabet_new_with_rank_func( const size_t size, const char *letters,
                                       char_rank_func rfunc );


/**
 * @brief Creates an alphabet from a string with letters in lexicographic order.
 * The ranks of letters will be put on a map so that they can be retrieved
 * in constant time.
 * @param size Number of letters.
 * @param letters String with letters in lexicographic order.
 */
alphabet *alphabet_new_with_rank_map(const size_t size, const char *letters);


/**
 * @brief Destructor.
 */
void alphabet_free(alphabet *ab);


/**
 * @brief Prints an alphabet representation to the standard output.
 */
void ab_print(alphabet *ab);


/**
 * @brief Returns the number of letters.
 */
size_t ab_size(alphabet *ab);


/**
 * @brief Indicates whether alphabet @p ab contains the character @c.
 */
bool ab_contains(alphabet *ab, char c); 

/**
 * @brief Returns the letter of given rank. If rank >= alphabet size, the
 * behaviour is undefined.
 */
char ab_char(alphabet *ab, size_t rank);


/**
 * @brief Returns the rank of a given character c.
 */
size_t ab_rank(alphabet *ab, char c);


/**
 * @brief Returns individual letter counts of a given string.
 */
size_t *ab_count(alphabet *ab, char *str, size_t slen);


/**
 * @brief Compute cumulative char frequencies
 * cumul_count(T, c) = Sum for a<c count(T,a)
 *
 * Example:
 *
 *        T =abracadabra   AB={a,b,c,d,r}
 *
 *           i    :   0    1    2    3    4    5
 *     alphabet   :   a    b    c    d    r
 *      chr_cnt   :   5    2    1    1    2
 *  cumul_freqs   :   0    5    7    8    9   11
 *
 *
 * Therefore, if the alphabet is A=<a[0], ..., a[s-1]>, then
 * the number of chars in the interval a[i]...a[j-1] is given by
 * cumul_count[j] - cumul_count[i]
 */
size_t  *ab_cumul_count(alphabet *ab, char *str, size_t slen);


/**
 * @brief Returns individual letter counts of a given stream.
 * @see ab_count
 */
size_t *ab_count_stream(alphabet *ab, strstream *sst);


/**
 * @brief Compute cumulative char frequencies of a given stream
 */
size_t  *ab_cumul_count_stream(alphabet *ab, strstream *sst);


/**
 * @brief Return the string representation of an alphabet, i.e. a string
 * composed of its symbols in lexicographic order.
 */
const char *ab_as_str(alphabet *ab);



///**
// * @brief Marshall an alphabet object
// */
//void ab_marshall(marshallctx *ctx, alphabet *ab);
//
//
///**
// * @brief Unmarshall and return an alphabet object
// */
//alphabet *ab_unmarshall(marshallctx *ctx);


#endif
