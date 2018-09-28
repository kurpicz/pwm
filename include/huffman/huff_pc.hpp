/*******************************************************************************
 * include/huffman/huff_pc.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

template <typename AlphabetType, typename ContextType, typename HuffCodes>
void huff_pc(AlphabetType const* text,
             uint64_t const size,
             uint64_t const levels,
             HuffCodes const& codes,
             ContextType& ctx)
{
  auto& zeros = ctx.zeros();
  auto& borders = ctx.borders();
  auto& bv = ctx.bv();

  // While calculating the histogram, we also compute the first level
  {
    ctx.hist(0, 0) = size;

    auto process_symbol = [&](auto& word, auto pos) {
      const code_pair cp = codes.encode_symbol(text[pos]);
      word <<= 1;
      word |= cp[0];
      for (size_t level = 1; level <= cp.code_length; level++) {
        auto prefix = cp.prefix(level);
        ctx.hist(level, prefix)++;
      }
    };

    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (uint64_t i = 0; i < 64; ++i) {
        process_symbol(word, cur_pos + i);
      }
      bv[0][cur_pos >> 6] = word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (uint64_t i = 0; i < size - cur_pos; ++i) {
        process_symbol(word, cur_pos + i);
      }
      word <<= (64 - (size & 63ULL));
      bv[0][size >> 6] = word;
    }
  }

  // Now we compute the WX top-down, since the histograms are already computed
  for (uint64_t level = 1; level < levels; level++) {
    uint64_t blocks = 1ull << level;

    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    borders[0] = 0;
    for (uint64_t i = 1; i < blocks; ++i) {
      auto const prev_block = ctx.rho(level, i - 1);
      auto const this_block = ctx.rho(level, i);

      borders[this_block] = borders[prev_block] + ctx.hist(level, prev_block);
      // NB: The above calulcation produces _wrong_ border offsets
      // for huffman codes that are one-shorter than the current level.
      //
      // Since those codes will not be used in the loop below, this does not
      // produce wrong or out-of-bound accesses.

      DCHECK(!ContextType::compute_rho);
      // TODO: Due to the forward iteration, this can not currently work
      /*
      if (ContextType::compute_rho)  {
        ctx.set_rho(level - 1, i - 1, prev_block >> 1);
      }
      */
    }

    if (ContextType::compute_zeros) {
      // If we compute zeros, we are working on a WM instead of a WT.
      // For a WM, borders is permuted with rho such that
      // borders[1] contains the position of the first 1-bit block.
      zeros[level - 1] = borders[1];
    }

    // Now we insert the bits with respect to their bit prefixes
    for (uint64_t i = 0; i < size; ++i) {
      const code_pair cp = codes.encode_symbol(text[i]);
      if (level < cp.code_length) {
        uint64_t const prefix = cp.prefix(level);
        uint64_t const pos = borders[prefix]++;
        uint64_t const bit = cp[level];
        uint64_t const word_pos = 63ULL - (pos & 63ULL);
        bv[level][pos >> 6] |= (bit << word_pos);
      }
    }
  }
}

/******************************************************************************/
