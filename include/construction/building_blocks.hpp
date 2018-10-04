/*******************************************************************************
 * include/huffman/huff_pc.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

template<typename text_t, typename ctx_t, typename bv_t>
inline void scan_text_compute_first_level_bv_and_last_level_hist(
  text_t const& text,
  size_t const size,
  uint64_t const levels,
  bv_t& bv,
  ctx_t& ctx
) {
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
}
