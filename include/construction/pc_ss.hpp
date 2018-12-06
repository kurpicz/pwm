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
void pc_ss(AlphabetType const* const text,
           uint64_t const size,
           uint64_t const levels,
           ContextType& ctx) {

  auto& bv = ctx.bv();

  scan_text_compute_first_level_bv_and_last_level_hist(text, size, levels, bv,
                                                       ctx);

  bottom_up_compute_hist_borders_optional_zeros_rho(size, levels, ctx);

  for (uint64_t i = 0; i < size; ++i) {
    auto const c = text[i];
    for (uint64_t level = levels - 1; level > 0; --level) {
      auto&& borders = ctx.borders_at_level(level);
      write_symbol_bit(bv, level, levels, borders, c);
    }
  }
}

/******************************************************************************/
