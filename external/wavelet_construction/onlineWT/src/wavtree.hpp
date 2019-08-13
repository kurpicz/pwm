#ifndef  ONLINE_WAVTREE_H
#define  ONLINE_WAVTREE_H

#include <stddef.h>

#include "alphabet.hpp"

/**
 * @file onlinewavelettree.h
 * @author Paulo Fonseca, Israel Batista
 *
 * @brief Wavelet tree ADT.
 *
 * The Wavelet Tree is a succinct data structure to store strings in
 * compressed space. It generalizes the rank and select
 * operations defined on bitvectors to arbitrary alphabets.
 * (https://en.wikipedia.org/wiki/Wavelet_Tree)
 */


/**
 * @brief Wavelet tree shape
 */
typedef enum {
    WT_BALANCED = 0, /**< Balanced, i.e. at each node the alphabet is split in halves. */
    WT_HUFFMAN = 1   /**< Shaped after a Huffman code tree */
}
wtshape;

typedef struct _wavtree wavtree;


/**
 * @brief Creates a wavelet tree representing a string @p str over an alphabet
 *        @p ab.
 * @param ab The base alphabet.
 * @param str The indexed string.
 * @param len The length of the indexed string.
 * @param mode The WT layout.
 */
wavtree *wavtree_new(alphabet *ab, char *str, size_t len, wtshape shape);


/**
 * @brief Creates a wavelet tree from a stream with known alphabet @p ab.
 * @param ab The base alphabet.
 * @param sst The source stream.
 * @param mode The WT layout.
 */
wavtree *wavtree_new_from_stream(alphabet *ab, strstream *sst, wtshape shape);


/**
 * @brief Create a wavelet tree from a stream with unknown alphabet.
 * @param sst The source stream.
 * @param mode The WT layout.
 */
wavtree *wavtree_new_online_from_stream(strstream *sst, wtshape shape);


/**
 * @brief Destructor
 */
void wavtree_free(wavtree *wt);


/**
 * @brief Computes the rank of a given position @p pos, defined as
 *        rank[pos] = # of positions 0<=j<=pos s.t. str[j]==str[pos], where
 *        str is string represented by the WT.
 */
size_t wavtree_rank_pos(wavtree *wt, size_t pos);


/**
 * @brief Computes the @pos c-rank of a given position @p pos, defined as
 *        the # of positions 0<=j<=@p pos s.t. str[j]==@p c, where
 *        str is string represented by the WT.
 */
size_t wavtree_rank(wavtree *wt, size_t pos, char c);


/**
 * @brief Computes the position of the occurence of a char @p c with
 *        a given @p rank in the string represented by the WT.
 *        select(c,r) = j s.t. str[j]==c and rank[j]==r
 */
size_t wavtree_select(wavtree *wt, char c, size_t rank);


/**
 * @brief Returns the char at position @p pos in the string represented by
 *        the WT. Notice that the WT does <b>not</b> explicitly store
 *        the string.
 */
char wavtree_char(wavtree *wt, size_t pos);


/**
 * @brief Prints a representation of the WT to standard output.
 */
void wavtree_print(wavtree *wt);

#endif
