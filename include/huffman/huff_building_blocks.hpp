/*******************************************************************************
 * include/huffman/huff_building_blocks.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "construction/building_blocks.hpp"

template <typename text_t, typename ctx_t, typename bv_t, typename codes_t>
inline void
huff_scan_text_compute_first_level_bv_and_full_hist(text_t const& text,
                                                    size_t const size,
                                                    bv_t& bv,
                                                    ctx_t& ctx,
                                                    codes_t const& codes) {
  ctx.hist(0, 0) = size;

  auto process_symbol = [&](auto& word, auto pos) {
    const code_pair cp = codes.encode_symbol(text[pos]);
    word <<= 1;
    word |= cp[0];
    // TODO: Try to find more efficient scheme that does not
    // do O(level) individual add operations.
    //
    // Eg, just adding at code_length, and summing all up afterwards?
    for (size_t level = 1; level <= cp.code_length; level++) {
      auto prefix = cp.prefix(level);
      ctx.hist(level, prefix)++;
    }
  };

  uint64_t cur_pos = 0;
  for (; cur_pos + 64 <= size; cur_pos += 64) {
    uint64_t word = 0ULL;
    for (uint64_t i = 0; i < 64; ++i) {
      process_symbol(word, cur_pos + i);
    }
    bv[0][cur_pos >> 6] = word;
  }
  if (size & 63ULL) {
    uint64_t word = 0ULL;
    for (uint64_t i = 0; i < size - cur_pos; ++i) {
      process_symbol(word, cur_pos + i);
    }
    word <<= (64 - (size & 63ULL));
    bv[0][size >> 6] = word;
  }
}
