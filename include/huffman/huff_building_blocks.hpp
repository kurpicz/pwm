/*******************************************************************************
 * include/huffman/huff_building_blocks.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 * Copyright (C) 2019 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <algorithm>
#include <unordered_map>

#include "construction/building_blocks.hpp"
#include "util/bit_reverse.hpp"

template <typename text_t, typename ctx_t, typename bv_t, typename codes_t>
inline void
huff_scan_text_compute_first_level_bv_and_full_hist(text_t const& text,
                                                    size_t const size,
                                                    bv_t& bv,
                                                    ctx_t& ctx,
                                                    codes_t const& codes) {
  ctx.hist_at_level(0)[0] = size;

  auto process_symbol = [&](auto& word, auto pos) {
    const code_pair cp = codes.encode_symbol(text[pos]);
    word <<= 1;
    word |= cp[0];
    // TODO: Try to find more efficient scheme that does not
    // do O(level) individual add operations.
    //
    // Eg, just adding at code_length, and summing all up afterwards?
    for (size_t level = 1; level <= cp.code_length; level++) {
      auto&& hist = ctx.hist_at_level(level);
      auto prefix = cp.prefix(level);
      hist[prefix]++;
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

template <typename ctx_t, typename borders_t>
inline void huff_compute_borders_optional_zeros_rho(uint64_t level,
                                                    uint64_t blocks,
                                                    ctx_t& ctx,
                                                    borders_t&& borders) {
  auto&& hist = ctx.hist_at_level(level);

  // Compute the starting positions of characters with respect to their
  // bit prefixes and the bit-reversal permutation
  borders[0] = 0;
  uint64_t containing_prev_block = 0;
  for (uint64_t i = 1; i < blocks; ++i) {
    // This is only required if we compute matrices
    [[maybe_unused]] auto const prev_block = ctx.rho(level, i - 1);
    auto const this_block = ctx.rho(level, i);

    if (hist.count(this_block) > 0) {
      borders[this_block] = borders[containing_prev_block] +
        hist[containing_prev_block];
      containing_prev_block = this_block;
    }
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

template <typename ctx_t, typename borders_t>
inline void huff_compute_borders_optional_zeros_rho(uint64_t level,
                                                    ctx_t& ctx,
                                                    borders_t&& borders) {
  auto&& hist = ctx.hist_at_level(level);
  std::vector<uint64_t> keys;
  keys.reserve(hist.size());
  std::for_each(hist.begin(), hist.end(), [&](auto const item) {
    keys.push_back(item.first);
  });

  if constexpr (ctx_t::compute_zeros) /* is wavelet matrix */ {
    std::sort(keys.begin(), keys.end(), [](uint64_t const a, uint64_t const b) {
      return reverse_bits(a) < reverse_bits(b);
    });
  } else /* is wavelet tree */{
    std::sort(keys.begin(), keys.end());
  }

  // Compute the starting positions of characters with respect to their
  // bit prefixes and the bit-reversal permutation
  std::cout << "keys.size() " << keys.size() << std::endl;
  borders.clear();
  borders[0] = 0;
  for (uint64_t i = 1; i < keys.size(); ++i) {
    borders[i] = borders[i - 1] + hist[i - 1];
  }

  if constexpr (ctx_t::compute_zeros) {
    // If we compute zeros, we are working on a WM instead of a WT.
    // For a WM, borders is permuted with rho such that
    // borders[1] contains the position of the first 1-bit block.
    ctx.zeros()[level - 1] = borders[1];
  }
}


/******************************************************************************/
