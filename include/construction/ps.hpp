/*******************************************************************************
 * include/util/ps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "construction/building_blocks.hpp"

template <typename AlphabetType, typename ContextType>
void ps(AlphabetType const* const text,
        uint64_t const size,
        uint64_t const levels,
        ContextType& ctx) {
  auto sorted_text_ = std::vector<AlphabetType>(size);
  auto sorted_text = span<AlphabetType>(sorted_text_);
  uint64_t cur_alphabet_size = (1 << levels);

  auto&& zeros = ctx.zeros();
  auto& bv = ctx.bv();

  scan_text_compute_first_level_bv_and_last_level_hist(text, size, levels, bv,
                                                       ctx);

  // The number of 0s at the last level is the number of "even" characters
  if (ContextType::compute_zeros) {
    auto&& hist = ctx.hist_at_level(levels);
    for (uint64_t i = 0; i < cur_alphabet_size; i += 2) {
      zeros[levels - 1] += hist[i];
    }
  }

  // Now we compute the WM bottom-up, i.e., the last level first
  for (uint64_t level = levels - 1; level > 0; --level) {
    auto&& borders = ctx.borders_at_level(level);
    auto&& hist = ctx.hist_at_level(level);
    auto&& next_hist = ctx.hist_at_level(level + 1);

    // Update the maximum value of a feasible a bit prefix and update the
    // histogram of the bit prefixes
    cur_alphabet_size >>= 1;
    for (uint64_t i = 0; i < cur_alphabet_size; ++i) {
      hist[i] = next_hist[i << 1] + next_hist[(i << 1) + 1];
    }

    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    compute_borders_optional_zeros_rho(level, cur_alphabet_size, ctx, borders);

    // Now we sort the text utilizing counting sort and the starting positions
    // that we have computed before
    for (uint64_t i = 0; i < size; ++i) {
      const auto cur_char = text[i];
      sorted_text[borders[cur_char >> (levels - level)]++] = cur_char;
    }

    // Since we have sorted the text, we can simply scan it from left to right
    // and for the character at position $i$ we set the $i$-th bit in the
    // bit vector accordingly
    write_bits_wordwise(0, size, bv[level], [&](uint64_t const i) {
      uint64_t const bit = ((sorted_text[i] >> ((levels - 1) - level)) & 1ULL);
      return bit;
    });
  }

  ctx.hist_at_level(0)[0] = size;
}

/******************************************************************************/
