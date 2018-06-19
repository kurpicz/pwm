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

template <typename AlphabteType, bool is_tree_, bool is_semi_external = false>
class wx_pc_ss {

public:
  static constexpr bool is_parallel = false;
  static constexpr bool is_tree     = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabteType);
  static constexpr bool  is_huffman_shaped = false;

  using ctx_t = ctx_compute_borders<is_tree, is_semi_external>;

  template <typename InputType>
  static wavelet_structure<is_semi_external> compute(const InputType& text, const uint64_t size,
    const uint64_t levels) {

    if (size == 0) {
      return wavelet_structure<is_semi_external>();
    }

    const auto rho = rho_dispatch<is_tree>::create(levels);
    auto ctx = ctx_t(size, levels, rho);
    
    pc_ss(text, size, levels, ctx);

    if (ctx_t::compute_zeros) {
      return wavelet_structure<is_semi_external>(std::move(ctx.bv()), std::move(ctx.zeros()));
    } else {
      return wavelet_structure<is_semi_external>(std::move(ctx.bv()));
    }
  }
}; // class wc_pc_ss

/******************************************************************************/
