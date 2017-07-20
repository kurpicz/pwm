/*******************************************************************************
 * include/prefix_counting.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef PREFIX_COUNTING_HEADER
#define PREFIX_COUNTING_HEADER

#include <cstring>
#include <vector>

#include "common.hpp"
#include "pc.hpp"

namespace construction {

template <typename AlphabetType,
      typename SizeType,
      template <typename> class WaveletStructure>
static WaveletStructure<SizeType> prefix_counting(const AlphabetType* text,
  const SizeType size, const SizeType levels) {
  WaveletStructure<SizeType> result(levels);

  if (size == 0) {
    return result;
  }
  SizeType cur_max_char = (1 << levels);
  std::vector<SizeType> perm = WaveletStructure<SizeType>::permutation(levels - 1);
  std::vector<SizeType> borders(cur_max_char, 0);

  auto ctx = SingleThreaded<SizeType> {
    levels, cur_max_char
  };

  result.bit_vector[0] = new uint64_t[word_size(size)];
  // memset is ok (all to 0)
  memset(result.bit_vector[0], 0, (word_size(size)) * sizeof(uint64_t));
  // While initializing the histogram, we also compute the first level
  SizeType cur_pos = 0;
  for (; cur_pos + 64 <= size; cur_pos += 64) {
    uint64_t word = 0ULL;
    for (SizeType i = 0; i < 64; ++i) {
      ++ctx.hist(levels, text[cur_pos + i]);
      word <<= 1;
      word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
    }
    result.bit_vector[0][cur_pos >> 6] = word;
  }
  if (size & 63ULL) {
    uint64_t word = 0ULL;
    for (SizeType i = 0; i < size - cur_pos; ++i) {
      ++ctx.hist(levels, text[cur_pos + i]);
      word <<= 1;
      word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
    }
    word <<= (64 - (size & 63ULL));
    result.bit_vector[0][size >> 6] = word;
  }

  // The number of 0s at the last level is the number of "even" characters
  if (WaveletStructure<SizeType>::is_tree) {
    for (SizeType i = 0; i < cur_max_char; i += 2) {
      result.zeros[levels - 1] += ctx.hist(levels, i);
    }
  }

  // Now we compute the WM bottom-up, i.e., the last level first
  for (SizeType level = levels - 1; level > 0; --level) {
    const SizeType prefix_shift = (levels - level);
    const SizeType cur_bit_shift = prefix_shift - 1;

    result.bit_vector[level] = new uint64_t[word_size(size)];
    // memset is ok (all to 0)
    memset(result.bit_vector[level], 0, (word_size(size)) * sizeof(uint64_t));

    // Update the maximum value of a feasible a bit prefix and update the
    // histogram of the bit prefixes
    cur_max_char >>= 1;
    for (SizeType i = 0; i < cur_max_char; ++i) {
      ctx.hist(level, i)
        = ctx.hist(level + 1, i << 1) + ctx.hist(level + 1, (i << 1) + 1);
    }

    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    borders[0] = 0;
    for (SizeType i = 1; i < cur_max_char; ++i) {
      borders[perm[i]] = borders[perm[i - 1]] +
      ctx.hist(level, perm[i - 1]);
      perm[i - 1] >>= 1;
    }
    // The number of 0s is the position of the first 1 in the previous level
    if (WaveletStructure<SizeType>::is_tree) {
      result.zeros[level - 1] = borders[1];
    }

    // Now we insert the bits with respect to their bit prefixes
    for (SizeType i = 0; i < size; ++i) {
      const SizeType pos = borders[text[i] >> prefix_shift]++;
      result.bit_vector[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
      << (63ULL - (pos & 63ULL)));
    }
  }

  return result;
}

} // namespace construction

#endif // PREFIX_COUNTING_HEADER

/******************************************************************************/
