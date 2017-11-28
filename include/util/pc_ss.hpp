/*******************************************************************************
 * include/util/pc_ss.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

template <typename InputType, typename ContextType>
void pc_ss(const InputType& text, uint64_t const size, const uint64_t levels,
  ContextType& ctx) {

  auto& bv = ctx.bv().raw_data();

  // While initializing the histogram, we also compute the first level
  uint64_t cur_pos = 0;
  for (; cur_pos + 64 <= size; cur_pos += 64) {
    uint64_t word = 0ULL;
    for (uint64_t i = 0; i < 64; ++i) {
      ++ctx.hist(levels, text[cur_pos + i]);
      word <<= 1;
      word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
    }
    bv[0][cur_pos >> 6] = word;
  }
  if (size & 63ULL) {
    uint64_t word = 0ULL;
    for (uint64_t i = 0; i < size - cur_pos; ++i) {
      ++ctx.hist(levels, text[cur_pos + i]);
      word <<= 1;
      word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
    }
    word <<= (64 - (size & 63ULL));
    bv[0][size >> 6] = word;
  }

  ctx.fill_borders();

  for (uint64_t i = 0; i < size; ++i) {
    for (uint64_t level = levels - 1; level > 0; --level) {
      const uint64_t prefix_shift = (levels - level);
      const uint64_t cur_bit_shift = prefix_shift - 1;
      const uint64_t pos = ctx.borders(level, text[i] >> prefix_shift)++;
      bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
        << (63ULL - (pos & 63ULL)));
    }
  }

  if (levels > 1) { // TODO check condition
    ctx.hist(0, 0) = ctx.hist(1, 0) + ctx.hist(1, 1);
  }
}

/******************************************************************************/
