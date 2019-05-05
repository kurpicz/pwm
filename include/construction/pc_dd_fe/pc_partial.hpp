/*******************************************************************************
 * include/construction/pc_dd_fe/pc_partial.hpp
 *
 * Copyright (C) 2019 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "construction/building_blocks.hpp"
#include "construction/pc_dd_fe/ctx_partial.hpp"

// compute local wavelet tree with word aligned intervals
template <typename AlphabetType>
void pc_partial(AlphabetType const* text, ctx_partial * ctx) {

  const uint64_t size = ctx->size();
  const uint64_t levels = ctx->levels();

  uint64_t cur_alphabet_size = (1ULL << levels);

  auto& borders = ctx->borders();
  auto& hist = ctx->hist();
  auto& bv = ctx->bv();
  auto& data_size = ctx->data_size();
  auto& level_data_sizes = ctx->level_data_sizes();

  scan_text_compute_first_level_bv_and_last_level_hist(text, size, levels, bv,
                                                       *ctx);

  // Now we compute the WM bottom-up, i.e., the last level first
  for (uint64_t level = levels - 1; level > 0; --level) {
    const uint64_t cur_prefix_shift = (levels - level);
    const uint64_t cur_bit_shift = cur_prefix_shift - 1;

    const auto& prev_hist = hist[level + 1];
    const auto& cur_hist = hist[level];

    // Update the maximum value of a feasible bit prefix and update the
    // histogram of the bit prefixes
    cur_alphabet_size >>= 1;
    borders[0] = 0;
    cur_hist[0] = prev_hist[0] + prev_hist[1];

    for (uint64_t i = 1; i < cur_alphabet_size; ++i) {
      // borders are word aligned
      borders[i] = borders[i - 1] + (((cur_hist[i - 1] + 63) >> 6) << 6);
      cur_hist[i] = prev_hist[i << 1] + prev_hist[(i << 1) + 1];
    }
    level_data_sizes[level] = (borders[cur_alphabet_size - 1] >> 6) +
                              ((cur_hist[cur_alphabet_size - 1] + 63) >> 6);
    data_size += level_data_sizes[level];

    // Now we insert the bits with respect to their bit prefixes
    for (uint64_t i = 0; i < size; ++i) {
      const uint64_t pos = borders[text[i] >> cur_prefix_shift]++;
      bv[level][pos >> 6] |=
          (((text[i] >> cur_bit_shift) & 1ULL) << (63ULL - (pos & 63ULL)));
    }
  }

  level_data_sizes[0] = ((size + 63) >> 6);
  data_size += level_data_sizes[0];
  hist[0][0] = size;
}

/******************************************************************************/
