/*******************************************************************************
 * include/wx_pc_ss.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/ctx_compute_borders.hpp"
#include "util/wavelet_structure.hpp"

template <typename AlphabteType, bool is_matrix>
class wx_pc_ss {
  using ctx_t = ctx_compute_borders<is_matrix>;

public:
  static constexpr bool is_parallel = false;
  static constexpr bool is_tree     = !is_matrix;
  static constexpr uint8_t word_width = sizeof(AlphabteType);

  static wavelet_structure compute(AlphabteType const* const text,
                                   const uint64_t size,
                                   const uint64_t levels) {
    if (size == 0) {
      return wavelet_structure();
    }

    const auto rho = rho_dispatch<is_matrix>::create(levels);
    auto ctx = ctx_t(size, levels, rho);
    auto& bv = ctx.bv().vec();

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

    if (ctx_t::compute_zeros) {
      return wavelet_structure(std::move(ctx.bv()), std::move(ctx.zeros()));
    } else {
      return wavelet_structure(std::move(ctx.bv()));
    }
  }
}; // class wc_pc_ss
