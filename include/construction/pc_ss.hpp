/*******************************************************************************
 * include/util/pc_ss.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "construction/building_blocks.hpp"

template <typename AlphabetType, typename ContextType>
void pc_ss(AlphabetType const* const text, uint64_t const size,
  uint64_t const levels, ContextType& ctx) {

  auto& bv = ctx.bv();

  scan_text_compute_first_level_bv_and_last_level_hist(
    text, size, levels, bv, ctx
  );

  bottom_up_compute_hist_and_borders_and_optional_zeros(levels, ctx);

  for (uint64_t i = 0; i < size; ++i) {
    for (uint64_t level = levels - 1; level > 0; --level) {
      const uint64_t prefix_shift = (levels - level);
      const uint64_t cur_bit_shift = prefix_shift - 1;
      const uint64_t pos = ctx.borders(level, text[i] >> prefix_shift)++;
      bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
        << (63ULL - (pos & 63ULL)));
    }
  }
}

/******************************************************************************/
