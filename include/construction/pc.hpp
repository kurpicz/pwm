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
  uint64_t cur_alphabet_size = (1 << levels);

  auto&& zeros = ctx.zeros();
  auto& bv = ctx.bv();

  scan_text_compute_first_level_bv_and_last_level_hist(text, size, levels, bv,
                                                       ctx);

  // The number of 0s at the last level is the number of "even" characters
  if (ContextType::compute_zeros) {
    auto hist = ctx.hist_at_level(levels);
    for (uint64_t i = 0; i < cur_alphabet_size; i += 2) {
      zeros[levels - 1] += hist[i];
    }
  }

  // Now we compute the WM bottom-up, i.e., the last level first
  for (uint64_t level = levels - 1; level > 0; --level) {
    const uint64_t prefix_shift = (levels - level);
    const uint64_t cur_bit_shift = prefix_shift - 1;
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
    compute_borders_and_optional_zeros_and_optional_rho(level,
                                                        cur_alphabet_size, ctx);
    auto&& borders = ctx.borders_at_level(level);

    // Now we insert the bits with respect to their bit prefixes
    for (uint64_t i = 0; i < size; ++i) {
      const uint64_t pos = borders[text[i] >> prefix_shift]++;
      bv[level][pos >> 6] |=
          (((text[i] >> cur_bit_shift) & 1ULL) << (63ULL - (pos & 63ULL)));
    }
  }

  ctx.hist_at_level(0)[0] = size;
}

/******************************************************************************/
