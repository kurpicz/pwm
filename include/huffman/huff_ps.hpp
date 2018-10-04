/*******************************************************************************
 * include/huffman/huff_ps.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstdint>
#include <cstddef>
#include <iostream>

#include "huffman/huff_codes.hpp"
#include "huffman/huff_building_blocks.hpp"

template <typename AlphabetType, typename ContextType, typename HuffCodes>
void huff_ps(AlphabetType const* text,
             uint64_t const size,
             uint64_t const levels,
             HuffCodes const& codes,
             ContextType& ctx,
             AlphabetType* const sorted_text,
             std::vector<uint64_t> const& level_sizes)
{
  auto& borders = ctx.borders();
  auto& bv = ctx.bv();

  // While calculating the histogram, we also compute the first level
  huff_scan_text_compute_first_level_bv_and_full_hist(
    text, size, bv, ctx, codes
  );

  // Now we compute the WX top-down, since the histograms are already computed
  for (uint64_t level = levels - 1; level > 0; --level) {
    uint64_t blocks = 1ull << level;

    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    compute_borders_and_optional_zeros(level, blocks, ctx);

    // Now we sort the text utilizing counting sort and the starting positions
    // that we have computed before
    for (uint64_t i = 0; i < size; ++i) {
      auto const cur_char = text[i];
      const code_pair cp = codes.encode_symbol(cur_char);

      // TODO: Make use of previously reduced sorted_text to
      // reduce iteration time?
      if (level < cp.code_length) {
        uint64_t const prefix = cp.prefix(level);
        uint64_t const pos = borders[prefix]++;
        sorted_text[pos] = cur_char;
      }
    }

    uint64_t const level_size = level_sizes[level];

    // Now we insert the bits with respect to their bit prefixes
    write_bits_wordwise(0, level_size, bv[level], [&](uint64_t const i) {
      const code_pair cp = codes.encode_symbol(sorted_text[i]);
      DCHECK(level < cp.code_length);
      uint64_t const bit = cp[level];
      return bit;
    });
  }
}

/******************************************************************************/
