/*******************************************************************************
 * include/util/ps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "construction/building_blocks.hpp"

template <typename AlphabetType, typename ContextType>
void ps(AlphabetType const* const text, uint64_t const size,
  uint64_t const levels, ContextType& ctx, AlphabetType* const sorted_text) {
  uint64_t cur_alphabet_size = (1 << levels);

  auto& zeros = ctx.zeros();
  auto& borders = ctx.borders();
  auto& bv = ctx.bv();

  scan_text_compute_first_level_bv_and_last_level_hist(
    text, size, levels, bv, ctx
  );

  // The number of 0s at the last level is the number of "even" characters
  if (ContextType::compute_zeros) {
    for (uint64_t i = 0; i < cur_alphabet_size; i += 2) {
      zeros[levels - 1] += ctx.hist(levels, i);
    }
  }

  // Now we compute the WM bottom-up, i.e., the last level first
  for (uint64_t level = levels - 1; level > 0; --level) {
    // Update the maximum value of a feasible a bit prefix and update the
    // histogram of the bit prefixes
    cur_alphabet_size >>= 1;
    for (uint64_t i = 0; i < cur_alphabet_size; ++i) {
      ctx.hist(level, i)
        = ctx.hist(level + 1, i << 1)
        + ctx.hist(level + 1, (i << 1) + 1);
    }

    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    borders[0] = 0;
    for (uint64_t i = 1; i < cur_alphabet_size; ++i) {
      auto const prev_rho = ctx.rho(level, i - 1);

      borders[ctx.rho(level, i)] =
        borders[prev_rho] + ctx.hist(level, prev_rho);

      if (ContextType::compute_rho)  {
        ctx.set_rho(level - 1, i - 1, prev_rho >> 1);
      }
    }

    // The number of 0s is the position of the first 1 in the previous level
    if (ContextType::compute_zeros) {
      zeros[level - 1] = borders[1];
    }

    // Now we sort the text utilizing counting sort and the starting positions
    // that we have computed before
    for (uint64_t i = 0; i < size; ++i) {
      const auto cur_char = text[i];
      sorted_text[borders[cur_char >> (levels - level)]++] = cur_char;
    }

    // Since we have sorted the text, we can simply scan it from left to right
    // and for the character at position $i$ we set the $i$-th bit in the
    // bit vector accordingly
    uint64_t cur_pos = 0;
    for (cur_pos = 0; cur_pos + 63 < size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (uint64_t i = 0; i < 64; ++i) {
        word <<= 1;
        word |= ((sorted_text[cur_pos + i] >> ((levels - 1) - level)) & 1ULL);
      }
      bv[level][cur_pos >> 6] = word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (uint64_t i = 0; i < size - cur_pos; ++i) {
        word <<= 1;
        word |= ((sorted_text[cur_pos + i] >> ((levels - 1) - level)) & 1ULL);
      }
      word <<= (64 - (size & 63ULL));
      bv[level][size >> 6] = word;
    }
  }

  ctx.hist(0, 0) = size;
}

/******************************************************************************/
