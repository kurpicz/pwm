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
void huff_pc_ss(AlphabetType const* text,
                uint64_t const size,
                uint64_t const levels,
                HuffCodes const& codes,
                ContextType& ctx,
                span<uint64_t const> const) {
  auto& bv = ctx.bv();

  // While calculating the histogram, we also compute the first level
  huff_scan_text_compute_first_level_bv_and_full_hist(text, size, bv, ctx,
                                                      codes);

  for (uint64_t level = levels - 1; level > 0; --level) {
    auto&& borders = ctx.borders_at_level(level);
    auto&& hist = ctx.hist_at_level(level);

    borders[0] = 0;
    uint64_t containing_prev_block = 0;
    for (uint64_t pos = 1; pos < ctx.hist_size(level); ++pos) {
      // This is only required if we compute matrices
      [[maybe_unused]] auto const prev_block = ctx.rho(level, pos - 1);
      auto const this_block = ctx.rho(level, pos);

      if (hist.count(this_block) > 0) {
        borders[this_block] = borders[containing_prev_block] +
          hist[containing_prev_block];
        containing_prev_block = this_block;
      }
      // NB: The above calulcation produces _wrong_ border offsets
      // for huffman codes that are one-shorter than the current level.
      //
      // Since those codes will not be used in the loop below, this does not
      // produce wrong or out-of-bound accesses.

      if (ContextType::compute_rho) {
        ctx.set_rho(level - 1, pos - 1, prev_block >> 1);
      }
    }

    // The number of 0s is the position of the first 1 in the previous level
    if constexpr (ContextType::compute_zeros) {
      ctx.zeros()[level - 1] = borders[1];
    }
  }

  for (uint64_t i = 0; i < size; ++i) {
    const code_pair cp = codes.encode_symbol(text[i]);
    for (uint64_t level = cp.code_length() - 1; level > 0; --level) {
      auto&& borders = ctx.borders_at_level(level);
      uint64_t const prefix = cp.prefix(level);
      uint64_t const pos = borders[prefix]++;
      uint64_t const bit = cp[level];
      uint64_t const word_pos = 63ULL - (pos & 63ULL);
      bv[level][pos >> 6] |= (bit << word_pos);
    }
  }
}

/******************************************************************************/
