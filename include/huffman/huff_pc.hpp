/*******************************************************************************
 * include/huffman/huff_pc.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstddef>
#include <cstdint>
#include <iostream>

#include "huffman/huff_building_blocks.hpp"
#include "huffman/huff_codes.hpp"

template <typename AlphabetType, typename ContextType, typename HuffCodes>
void huff_pc(AlphabetType const* text,
             uint64_t const size,
             uint64_t const levels,
             HuffCodes const& codes,
             ContextType& ctx,
             span<uint64_t const> const) {
  auto& bv = ctx.bv();

  // While calculating the histogram, we also compute the first level
  huff_scan_text_compute_first_level_bv_and_full_hist(text, size, bv, ctx,
                                                      codes);

  // Now we compute the WX top-down, since the histograms are already computed
  for (uint64_t level = levels - 1; level > 0; --level) {
    auto&& borders = ctx.borders_at_level(level);
    uint64_t blocks = 1ull << level;

    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    compute_borders_and_optional_zeros_and_optional_rho(level,
                                                        blocks, ctx, borders);

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
