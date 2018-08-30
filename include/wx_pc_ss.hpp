/*******************************************************************************
 * include/wx_pc_ss.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "construction/ctx_compute_borders.hpp"
#include "construction/pc_ss.hpp"
#include "construction/wavelet_structure.hpp"

template <typename AlphabteType, bool is_tree_>
class wx_pc_ss {

public:
  static constexpr bool is_parallel = false;
  static constexpr bool is_tree     = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabteType);
  static constexpr bool  is_huffman_shaped = false;

  using ctx_t = ctx_compute_borders<is_tree>;

  template <typename InputType>
  static wavelet_structure compute(const InputType& text, const uint64_t size,
    const uint64_t levels) {

    if(size == 0) {
      if constexpr (ctx_t::compute_zeros) { return wavelet_structure_matrix(); }
      else { return wavelet_structure_tree(); }
    }

    const auto rho = rho_dispatch<is_tree>::create(levels);
    auto ctx = ctx_t(size, levels, rho);

    pc_ss(text, size, levels, ctx);

    if constexpr (ctx_t::compute_zeros) {
      return wavelet_structure_matrix(
        std::move(ctx.bv()), std::move(ctx.zeros()));
    } else { return wavelet_structure_tree(std::move(ctx.bv())); }
  }
}; // class wc_pc_ss

/******************************************************************************/
