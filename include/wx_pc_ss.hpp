/*******************************************************************************
 * include/wx_pc_ss.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/ctx_compute_borders.hpp"
#include "util/pc_ss.hpp"
#include "util/wavelet_structure.hpp"

template <typename AlphabteType, bool is_matrix>
class wx_pc_ss {
  using ctx_t = ctx_compute_borders<is_matrix>;

public:
  static constexpr bool is_parallel = false;
  static constexpr bool is_tree     = !is_matrix;
  static constexpr uint8_t word_width = sizeof(AlphabteType);
  static constexpr bool  is_huffman_shaped = false;

  static wavelet_structure compute(AlphabteType const* const text,
                                   const uint64_t size,
                                   const uint64_t levels) {
    if (size == 0) {
      return wavelet_structure();
    }

    const auto rho = rho_dispatch<is_matrix>::create(levels);
    auto ctx = ctx_t(size, levels, rho);
    
    pc_ss(text, size, levels, ctx);

    if (ctx_t::compute_zeros) {
      return wavelet_structure(std::move(ctx.bv()), std::move(ctx.zeros()));
    } else {
      return wavelet_structure(std::move(ctx.bv()));
    }
  }
}; // class wc_pc_ss

/******************************************************************************/
