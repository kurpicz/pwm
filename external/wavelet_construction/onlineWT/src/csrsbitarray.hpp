#ifndef CSRSBITARRAY_H
#define CSRSBITARRAY_H

/**
 * @file csrsbitarray.h
 * @author Paulo Fonseca
 *
 * @brief Combined sampling rank&select static bitarray.
 * (https://www.dcc.uchile.cl/~gnavarro/ps/sea12.1.pdf)
 */

typedef struct _csrsbitarray csrsbitarray;


/**
 * @brief Creates a new r&s bitarray from raw source bits.
 * @param ba The source bitarray
 * @param len The array length
 */
csrsbitarray *csrsbitarr_new(byte_t *ba, size_t len);


/**
 * @brief Destructor.
 * @param free_data Indicates whether the source raw bitarray should be freed.
 */
void csrsbitarr_free(csrsbitarray *ba, bool free_data);


/**
 * @brief Prints a representations of the bitarray to standard output.
 */
void csrsbitarr_print(csrsbitarray *ba, size_t bytes_per_row);


/**
 * @brief Returns the bitarray's data.
 */
const byte_t *csrsbitarr_data(csrsbitarray *ba);

/**
 * @brief Returns the length of the bitarray.
 */
size_t csrsbitarr_len(csrsbitarray *ba);


/**
 * @brief Returns the bit at a certain position @p pos.
 */
byte_t csrsbitarr_get(csrsbitarray *ba, size_t pos);


/**
 * @brief Same as csrsbitarr_rank(@p ba, @p pos, 0).
 * @see csrsbitarr_rank
 */
size_t csrsbitarr_rank0(csrsbitarray *ba, size_t pos);


/**
 * @brief Same as csrsbitarr_rank(@p ba, @p pos, 1).
 * @see csrsbitarr_rank
 */
size_t csrsbitarr_rank1(csrsbitarray *ba, size_t pos);


/**
 * @brief Computes rank_@p bit(@p ba, @p pos) = # positions j<=@p pos
 * s.t. @p ba[j]==@p bit, for 0 <= @p pos < @p ba.len. If @pos>= @ba.len
 * returns the total number of positions with value == @p bit.
 */
size_t csrsbitarr_rank(csrsbitarray *ba, size_t pos, byte_t bit);


/**
 * @brief Same as csrsbitarr_select(@p ba, @p rank, 0).
 * @see csrsbitarr_select
 */
size_t csrsbitarr_select0(csrsbitarray *ba, size_t rank);


/**
 * @brief Same as csrsbitarr_select(@p ba, @p rank, 1).
 * @see csrsbitarr_select
 */
size_t csrsbitarr_select1(csrsbitarray *ba, size_t rank);


/**
 * @brief Computes select_@p bit(@p ba, @p rank) = j s.t.
 * @p ba[j]==@p bit and rank_@p bit(@p ba, j)=@p rank.
 * If no such position exists, return @p ba.len.
 */
size_t csrsbitarr_select(csrsbitarray *ba, size_t rank, byte_t bit);


/**
 * @brief Same as csrsbitarr_pred(@p ba, @p pos, 0).
 * @see csrsbitarr_pred
 */
size_t csrsbitarr_pred0(csrsbitarray *ba, size_t pos);


/**
 * @brief Same as csrsbitarr_pred(@p ba, @p pos, 1).
 * @see csrsbitarr_pred
 */
size_t csrsbitarr_pred1(csrsbitarray *ba, size_t pos);


/**
 * @brief Returns the rightmost position whose value is @p bit,
 * strictly to the left of @p pos, i.e max{j<pos | @p ba[j]==1}.
 * If no such position exists, returns @p ba.len.
 */
size_t csrsbitarr_pred(csrsbitarray *ba, size_t pos, byte_t bit);


/**
 * @brief Same as csrsbitarr_succ(@p ba, @p pos, 0).
 * @see csrsbitarr_succ
 */
size_t csrsbitarr_succ0(csrsbitarray *ba, size_t pos);


/**
 * @brief Same as csrsbitarr_succ(@p ba, @p pos, 1).
 * @see csrsbitarr_succ
 */
size_t csrsbitarr_succ1(csrsbitarray *ba, size_t pos);


/**
 * @brief Returns the leftmost position whose value is @p bit,
 * strictly to the right of @p pos, i.e min{j>pos | @p ba[j]==@p bit}.
 * If no such position exists, returns @p ba.len.
 */
size_t csrsbitarr_succ(csrsbitarray *ba, size_t pos, byte_t bit);


#endif
