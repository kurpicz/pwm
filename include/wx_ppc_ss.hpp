/*******************************************************************************
 * include/wx_ppc_ss.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/ctx_compute_borders.hpp"
#include "util/pc_ss.hpp"
#include "util/wavelet_structure.hpp"
#include "util/memory_types.hpp"

template <typename AlphabteType, bool is_tree_, memory_mode mem_mode_>
class wx_ppc_ss {

public:
  static constexpr bool is_parallel = true;
  static constexpr bool is_tree     = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabteType);
  static constexpr bool  is_huffman_shaped = false;
  static constexpr memory_mode mem_mode = mem_mode_;

  template <typename InputType, typename OutputType>
  static wavelet_structure<OutputType> compute(const InputType& text, const uint64_t size,
    const uint64_t levels) {

    using ctx_t = ctx_compute_borders<OutputType, is_tree>;

    if (size == 0) {
      return wavelet_structure<OutputType>();
    }

    const auto rho = rho_dispatch<is_tree>::create(levels);
    auto ctx = ctx_t(size, levels, rho);
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

    #pragma omp parallel num_threads(levels)
    {
      for (uint64_t i = 0; i < size; ++i) {
        #pragma omp for nowait
        for (uint64_t level = levels - 1; level > 0; --level) {
          const uint64_t prefix_shift = (levels - level);
          const uint64_t cur_bit_shift = prefix_shift - 1;
          const uint64_t pos = ctx.borders(level, text[i] >> prefix_shift)++;
          bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
            << (63ULL - (pos & 63ULL)));
        }
      }
    }

    if (ctx_t::compute_zeros) {
      return wavelet_structure<OutputType>(std::move(ctx.bv()), std::move(ctx.zeros()));
    } else {
      return wavelet_structure<OutputType>(std::move(ctx.bv()));
    }
  }
}; // class wx_ppc_ss

/******************************************************************************/
