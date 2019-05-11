/*******************************************************************************
 * include/huffman/huff_ps.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 * Copyright (C) 2019 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
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
void huff_ps(AlphabetType const* text,
             uint64_t const size,
             uint64_t const levels,
             HuffCodes const& codes,
             ContextType& ctx,
             span<uint64_t const> const level_sizes) {
  auto sorted_text_ = std::vector<AlphabetType>(size);
  auto sorted_text = span<AlphabetType>(sorted_text_);
  auto& bv = ctx.bv();

  std::vector<AlphabetType> mutable_text(text, text + size);

  // While calculating the histogram, we also compute the first level
  huff_scan_text_compute_first_level_bv_and_full_hist(text, size, bv, ctx,
                                                      codes);
  uint64_t text_length = size;
  for (uint64_t level = 1;
       level < std::min(levels, ContextType::UPPER_LEVEL_SIZE) &&
         bv[level].size() > 0; ++level) {
    auto&& borders = ctx.upper_borders_at_level(level);
    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    huff_compute_borders_optional_zeros_rho(level, ctx, borders);

    // Now we sort the text utilizing counting sort and the starting positions
    // that we have computed before
    uint64_t new_text_length = 0;
    for (uint64_t i = 0; i < text_length; ++i) {
      auto const cur_char = mutable_text[i];
      const code_pair cp = codes.encode_symbol(cur_char);

      // TODO: Make use of previously reduced sorted_text to
      // reduce iteration time?
      if (level < cp.code_length()) {
        uint64_t const prefix = cp.prefix(level);
        uint64_t const pos = borders[prefix]++;
        mutable_text[new_text_length++] = cur_char;
        sorted_text[pos] = cur_char;
      }
    }
    text_length = new_text_length;
    
    uint64_t const level_size = level_sizes[level];
    // Now we insert the bits with respect to their bit prefixes
    write_bits_wordwise(0, level_size, bv[level], [&](uint64_t const i) {
      code_pair const cp = codes.encode_symbol(sorted_text[i]);
      DCHECK(level < cp.code_length());
      uint64_t const bit = cp[level];
      return bit;
    });
  }
  for (uint64_t level = std::min(levels, ContextType::UPPER_LEVEL_SIZE);
       level < levels && bv[level].size() > 0; ++level) {
    auto&& borders = ctx.lower_borders_at_level(level);
    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    huff_compute_borders_optional_zeros_rho(level, ctx, borders);

    // Now we sort the text utilizing counting sort and the starting positions
    // that we have computed before
    uint64_t new_text_length = 0;
    for (uint64_t i = 0; i < text_length; ++i) {
      auto const cur_char = mutable_text[i];
      const code_pair cp = codes.encode_symbol(cur_char);

      // TODO: Make use of previously reduced sorted_text to
      // reduce iteration time?
      if (level < cp.code_length()) {
        uint64_t const prefix = cp.prefix(level);
        uint64_t const pos = borders[prefix]++;
        mutable_text[new_text_length++] = cur_char;
        sorted_text[pos] = cur_char;
      }
    }
    text_length = new_text_length;

    uint64_t const level_size = level_sizes[level];

    // Now we insert the bits with respect to their bit prefixes
    write_bits_wordwise(0, level_size, bv[level], [&](uint64_t const i) {
      const code_pair cp = codes.encode_symbol(sorted_text[i]);
      DCHECK(level < cp.code_length());
      uint64_t const bit = cp[level];
      return bit;
    });
  }
}

/******************************************************************************/
