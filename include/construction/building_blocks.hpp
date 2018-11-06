/*******************************************************************************
 * include/construction/building_blocks.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

template <typename level_bv_t, typename loop_body_t>
inline void write_bits_wordwise(uint64_t start,
                                uint64_t size,
                                level_bv_t const& level_bv,
                                loop_body_t body) {
  uint64_t cur_pos = start;
  for (; cur_pos + 64 <= size; cur_pos += 64) {
    uint64_t word = 0ULL;
    for (uint64_t i = 0; i < 64; ++i) {
      uint64_t const bit = body(cur_pos + i);
      word <<= 1;
      word |= bit;
    }
    level_bv[cur_pos >> 6] = word;
  }
  if (size & 63ULL) {
    uint64_t word = 0ULL;
    for (uint64_t i = cur_pos; i < size; ++i) {
      uint64_t const bit = body(i);
      word <<= 1;
      word |= bit;
    }
    word <<= (64 - (size & 63ULL));
    level_bv[size >> 6] = word;
  }
}

template <typename text_t, typename ctx_t, typename bv_t>
inline void
scan_text_compute_first_level_bv_and_last_level_hist(text_t const& text,
                                                     size_t const size,
                                                     uint64_t const levels,
                                                     bv_t& bv,
                                                     ctx_t& ctx) {
  write_bits_wordwise(0, size, bv[0], [&](uint64_t i) {
    ++ctx.hist(levels, text[i]);
    uint64_t const bit = ((text[i] >> (levels - 1)) & 1ULL);
    return bit;
  });
}

template <typename ctx_t>
inline void bottom_up_compute_hist_and_borders_and_optional_zeros(
    uint64_t const size, uint64_t const levels, ctx_t& ctx) {
  for (uint64_t level = levels - 1; level > 0; --level) {
    for (uint64_t pos = 0; pos < ctx.hist_size(level); ++pos) {
      ctx.hist(level, pos) =
          ctx.hist(level + 1, pos << 1) + ctx.hist(level + 1, (pos << 1) + 1);
    }

    ctx.borders(level, 0) = 0;
    for (uint64_t pos = 1; pos < ctx.hist_size(level); ++pos) {
      auto const prev_rho = ctx.rho(level, pos - 1);

      ctx.borders(level, ctx.rho(level, pos)) =
          ctx.borders(level, prev_rho) + ctx.hist(level, prev_rho);
    }

    // The number of 0s is the position of the first 1 in the previous level
    if constexpr (ctx_t::compute_zeros) {
      ctx.zeros()[level - 1] = ctx.borders(level, 1);
    }
  }
  ctx.hist(0, 0) = size;
}

template <typename ctx_t>
inline void compute_borders_and_optional_zeros_and_optional_rho(uint64_t level,
                                                                uint64_t blocks,
                                                                ctx_t& ctx) {
  auto& borders = ctx.borders();

  // Compute the starting positions of characters with respect to their
  // bit prefixes and the bit-reversal permutation
  borders[0] = 0;
  for (uint64_t i = 1; i < blocks; ++i) {
    auto const prev_block = ctx.rho(level, i - 1);
    auto const this_block = ctx.rho(level, i);

    borders[this_block] = borders[prev_block] + ctx.hist(level, prev_block);
    // NB: The above calulcation produces _wrong_ border offsets
    // for huffman codes that are one-shorter than the current level.
    //
    // Since those codes will not be used in the loop below, this does not
    // produce wrong or out-of-bound accesses.

    if (ctx_t::compute_rho) {
      ctx.set_rho(level - 1, i - 1, prev_block >> 1);
    }
  }

  if (ctx_t::compute_zeros) {
    // If we compute zeros, we are working on a WM instead of a WT.
    // For a WM, borders is permuted with rho such that
    // borders[1] contains the position of the first 1-bit block.
    ctx.zeros()[level - 1] = borders[1];
  }
}
