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
#include <chrono>
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

  std::vector<uint8_t> mutable_text(text, text + size);

  // While calculating the histogram, we also compute the first level
  huff_scan_text_compute_first_level_bv_and_full_hist(text, size, bv, ctx,
                                                      codes);

  // Now we compute the WX top-down, since the histograms are already computed
  uint64_t text_length = size;
  for (uint64_t level = 1; level < levels; ++level) {
    auto&& borders = ctx.borders_at_level(level);

    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    
    auto begin_time = std::chrono::high_resolution_clock::now();

    huff_compute_borders_optional_zeros_rho(level, ctx, borders);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = end_time - begin_time;
    auto millis =
      std::chrono::duration_cast<std::chrono::milliseconds>(duration);

    std::cout << "huffman_borders=" << millis.count() << std::endl;

    // Now we insert the bits with respect to their bit prefixes
    
    begin_time = std::chrono::high_resolution_clock::now();
    uint64_t write_pos = 0;
    for (uint64_t i = 0; i < text_length; ++i) {
      const code_pair cp = codes.encode_symbol(mutable_text[i]);
      if (level + 1 < cp.code_length()) {
        mutable_text[write_pos++] = mutable_text[i];
      }
      uint64_t const prefix = cp.prefix(level);
      uint64_t const pos = borders[prefix]++;
      uint64_t const bit = cp[level];
      uint64_t const word_pos = 63ULL - (pos & 63ULL);
      bv[level][pos >> 6] |= (bit << word_pos);
    }
    text_length = write_pos;
    end_time = std::chrono::high_resolution_clock::now();
    duration = end_time - begin_time;
    millis = std::chrono::duration_cast<std::chrono::milliseconds>(duration);

    std::cout << "huffman_fill=" << millis.count() << std::endl;
  }
}

/******************************************************************************/
