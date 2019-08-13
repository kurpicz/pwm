#ifndef HUFFCODE_H
#define HUFFCODE_H

#include <stddef.h>

#include "alphabet.hpp"
#include "bitsandbytes.hpp"
#include "cocadautil.hpp"
#include "strstream.hpp"

/**
 * @file huffcode.h
 * @author Paulo Fonseca
 * 
 * @brief Huffman coding data structures and algorithms.
 * 
 * A Huffman coding (HC) is a fixed-to-var prefix coding scheme that 
 * unambiguously associates every char c of an alphabet to a binary 
 * codeword h(c) such that if c!=c' then neither h(c) is a proper 
 * prefix of h(c') nor the opposite.
 * 
 * The HC is theoretically optimal in the information theoretic sense, meaning
 * that the average number of bits needed to encode a string T of length n 
 * is n*H(T) where H(T) denotes the 0-order entropy of T.
 * 
 * A HC over an alphabet A of length L is associated with a strictly binary 
 * trie - The Huffman Tree (HT) - with L leaves, each corresponding to a 
 * different character of A. The codeword of a char c is given by the label 
 * of the path from the root of the HT to its corresponding leaf.
 */


/**
 * Huffman code type
 */
typedef struct _huffcode huffcode;


/**
 * Huffman Tree node type
 */
typedef struct _hufftnode hufftnode;


/**
 * @brief Binary code. In case of prefix codes, such as the HC, it can be 
 *        unambiguously decoded in one pass from left to right.
 */
typedef struct {
    byte_t *rawcode; /**< raw bits */
    size_t code_len; /**< code length */
} bincode;


/**
 * @brief Creates a HC for an alphabet with associated letter frequencies.
 * @param ab The base alphabet.
 * @param freqs Individual letter frequencies in lexycographic order.
 */
huffcode *huffcode_new(alphabet *ab, size_t freqs[]);


/**
 * @brief Creates a HC for an alphabet from a source string.
 * @param ab The base alphabet.
 * @param src Source string from which letter frequencies are to be estimated.
 */
huffcode *huffcode_new_from_str(alphabet *ab, char *src);


/**
 * @brief Creates a HC for an alphabet from a source stream.
 * @param ab The base alphabet.
 * @param src Source stream from which letter frequencies are to be estimated.
 */
huffcode *huffcode_new_from_stream(alphabet *ab, strstream *src);


/**
 * @brief Creates a HC from a source stream.
 * @param src Source stream from which the support alphabet and 
 *        character frequencies are obtained.
 */
huffcode *huffcode_new_online_from_stream(strstream *sst);


/**
 * @brief Destructor.
 */
void huffcode_free(huffcode *hcode);


/**
 * @brief Prints a string representation of the HC to the std output. 
 */
void huffcode_print(huffcode *hcode);


/**
 * @brief Encodes a string into a binary code.
 */
bincode huffcode_encode(huffcode *hcode, char *str, size_t len);


/**
 * @brief Decodes a binary code to a string.
 */
char *huffcode_decode(huffcode *hcode, byte_t *code, size_t code_len);


/**
 * @brief Returns the (root of the) Huffman tree corresponding to a given HC. 
 */
hufftnode *huffcode_tree(huffcode *hcode);


/**
 * @brief Returns the alphabet from a given HC.
 */
alphabet *huffcode_ab(huffcode *code);


/**
 * @brief Checks whether a given HT node is a leaf.
 */
bool hufftnode_is_leaf(hufftnode *node);


/**
 * @brief Returns the left child of a HT node.
 */
hufftnode *hufftnode_left(hufftnode *node);


/**
 * @brief Returns the right child of a HT node.
 */
hufftnode *hufftnode_right(hufftnode *node);


/**
 * @brief Returns the alphabet mask of a HT @p node. The alphabet mask of a 
 *        @p node is a binary mask M[0..l-1] s.t. for each letter of the 
 *        alphabet a[j], M[j] = 1 iff the leaf corresponding to letter a[j] 
 *        is in the subtree rooted at this @p node.
 */
byte_t *hufftnode_ab_mask(hufftnode *node);


/**
 * @brief If the given HT @p node is a leaf, returns the lexycographic rank
 *        of the corresponding letter in the represented alphabet. 
 *        If the given HT @p node is not a leaf, returns the alphabet size.
 */
size_t hufftnode_char_rank(hufftnode *node);

#endif
