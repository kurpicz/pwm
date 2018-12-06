/*******************************************************************************
 * include/construction/building_blocks.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <omp.h>

#include "arrays/span.hpp"

// NB: This formulation as an inline-always function
// has been roghly tested to produce the least amount of machine code if
// used from `write_bits_wordwise()`.
template <typename level_bv_t, typename loop_body_t>
inline __attribute__((always_inline)) void
write_bits_word(uint64_t const start,
                uint64_t const size,
                level_bv_t const& level_bv,
                loop_body_t body) {
  DCHECK(size <= 64);
  uint64_t word = 0ULL;
  for (uint64_t i = 0; i < size; ++i) {
    uint64_t const bit = body(start + i);
    word <<= 1;
    word |= bit;
  }

  // NB. This gets optimized out if size is a constant 64
  word <<= (64 - size);

  level_bv[start >> 6] = word;
}

template <typename level_bv_t, typename loop_body_t>
inline void write_bits_wordwise(uint64_t start,
                                uint64_t size,
                                level_bv_t const& level_bv,
                                loop_body_t body) {
  for (uint64_t cur_pos = start; cur_pos + 64 <= size; cur_pos += 64) {
    write_bits_word(cur_pos, 64, level_bv, body);
  }
  uint64_t const remainder = size & 63ULL;
  if (remainder) {
    write_bits_word(size - remainder, remainder, level_bv, body);
  }
}

template <typename level_bv_t, typename loop_body_t>
inline void omp_write_bits_wordwise(uint64_t start,
                                    uint64_t size,
                                    level_bv_t const& level_bv,
                                    loop_body_t body) {
  const auto omp_rank = omp_get_thread_num();
  const auto omp_size = omp_get_num_threads();

  #pragma omp for
  for (int64_t scur_pos = start; scur_pos <= (int64_t(size) - 64);
       scur_pos += 64) {
    DCHECK(scur_pos >= 0);
    write_bits_word(scur_pos, 64, level_bv, body);
  }
  uint64_t const remainder = size & 63ULL;
  if (remainder && ((omp_rank + 1) == omp_size)) {
    write_bits_word(size - remainder, remainder, level_bv, body);
  }
}

template <typename text_t, typename ctx_t, typename bv_t>
inline void
scan_text_compute_first_level_bv_and_last_level_hist(text_t const& text,
                                                     size_t const size,
                                                     uint64_t const levels,
                                                     bv_t& bv,
                                                     ctx_t& ctx) {
  auto&& hist = ctx.hist_at_level(levels);
  write_bits_wordwise(0, size, bv[0], [&](uint64_t i) {
    ++hist[text[i]];
    uint64_t const bit = ((text[i] >> (levels - 1)) & 1ULL);
    return bit;
  });
}

template <typename ctx_t, typename borders_t>
inline void compute_borders_and_optional_zeros_and_optional_rho(
    uint64_t level, uint64_t blocks, ctx_t& ctx, borders_t&& borders) {
  auto&& hist = ctx.hist_at_level(level);

  // Compute the starting positions of characters with respect to their
  // bit prefixes and the bit-reversal permutation
  borders[0] = 0;
  for (uint64_t i = 1; i < blocks; ++i) {
    auto const prev_block = ctx.rho(level, i - 1);
    auto const this_block = ctx.rho(level, i);

    borders[this_block] = borders[prev_block] + hist[prev_block];
    // NB: The above calulcation produces _wrong_ border offsets
    // for huffman codes that are one-shorter than the current level.
    //
    // Since those codes will not be used in the loop below, this does not
    // produce wrong or out-of-bound accesses.

    if constexpr (ctx_t::compute_rho) {
      ctx.set_rho(level - 1, i - 1, prev_block >> 1);
    }
  }

  if constexpr (ctx_t::compute_zeros) {
    // If we compute zeros, we are working on a WM instead of a WT.
    // For a WM, borders is permuted with rho such that
    // borders[1] contains the position of the first 1-bit block.
    ctx.zeros()[level - 1] = borders[1];
  }
}

template <typename ctx_t>
inline void bottom_up_compute_hist_and_borders_and_optional_zeros(
    uint64_t const size, uint64_t const levels, ctx_t& ctx) {
  for (uint64_t level = levels - 1; level > 0; --level) {
    auto const blocks = ctx.hist_size(level);
    auto&& hist = ctx.hist_at_level(level);
    auto&& next_hist = ctx.hist_at_level(level + 1);
    auto&& borders = ctx.borders_at_level(level);

    for (uint64_t pos = 0; pos < blocks; ++pos) {
      hist[pos] = next_hist[pos << 1] + next_hist[(pos << 1) + 1];
    }

    compute_borders_and_optional_zeros_and_optional_rho(level,
                                                        blocks,
                                                        ctx,
                                                        borders);
  }
  ctx.hist_at_level(0)[0] = size;
}

template <typename bv_t, typename borders_t, typename alphabet_type>
inline __attribute__((always_inline)) void
write_symbol_bit(bv_t& bv,
                 uint64_t level,
                 uint64_t levels,
                 borders_t&& borders,
                 alphabet_type c) {
  // NB: The computations for `prefix_shift` and `cur_bit_shift` will most
  // likely be hoisted out of a inner loop by a optimizing compiler
  const uint64_t prefix_shift = (levels - level);
  const uint64_t cur_bit_shift = prefix_shift - 1;
  const uint64_t pos = borders[c >> prefix_shift]++;
  const uint64_t bit =
      (((c >> cur_bit_shift) & 1ULL) << (63ULL - (pos & 63ULL)));
  bv[level][pos >> 6] |= bit;
}
