/*******************************************************************************
 * include/util/pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "construction/building_blocks.hpp"

template <typename AlphabetType, typename ContextType>
void pc(AlphabetType const* text,
        const uint64_t size,
        const uint64_t levels,
        ContextType& ctx) {
  uint64_t cur_alphabet_size = (1ull << levels);

  auto&& zeros = ctx.zeros();
  auto& bv = ctx.bv();

  scan_text_compute_first_level_bv_and_last_level_hist(text, size, levels, bv,
                                                       ctx);

  // The number of 0s at the last level is the number of "even" characters
  if constexpr (ContextType::compute_zeros) {
    auto hist = ctx.hist_at_level(levels);
    for (uint64_t i = 0; i < cur_alphabet_size; i += 2) {
      zeros[levels - 1] += hist[i];
    }
  }

  // Now we compute the WM bottom-up, i.e., the last level first
  for (uint64_t level = levels - 1; level > 0; --level) {
    auto&& hist = ctx.hist_at_level(level);
    auto&& next_hist = ctx.hist_at_level(level + 1);

    // Update the maximum value of a feasible a bit prefix and update the
    // histogram of the bit prefixes
    cur_alphabet_size >>= 1;
    for (uint64_t i = 0; i < cur_alphabet_size; ++i) {
      hist[i] = next_hist[i << 1] + next_hist[(i << 1) + 1];
    }

    auto&& borders = ctx.borders_at_level(level);

    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    compute_borders_and_optional_zeros_and_optional_rho(level,
                                                        cur_alphabet_size,
                                                        ctx, borders);

    // Now we insert the bits with respect to their bit prefixes
    for (uint64_t i = 0; i < size; ++i) {
      write_symbol_bit(bv, level, levels, borders, text[i]);
    }
  }

  ctx.hist_at_level(0)[0] = size;
}

/******************************************************************************/
