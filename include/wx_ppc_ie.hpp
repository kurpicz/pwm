/*******************************************************************************
 * include/wx_pc_ie.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>

#include "arrays/memory_types.hpp"
#include "construction/ctx_single_level.hpp"
#include "construction/pc_semi_external.hpp"
#include "construction/wavelet_structure.hpp"

#include "wx_base.hpp"

template <typename AlphabetType, bool is_tree_>
class wx_ppc_ie : public wx_in_out_external<true, false> {

public:
  static constexpr bool is_parallel = true;
  static constexpr bool is_tree = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = false;

  template <typename InputType>
  static wavelet_structure
  compute(const InputType& text, const uint64_t size, const uint64_t levels) {

    using ctx_t = ctx_single_level<is_tree>;

    if (size == 0) {
      if constexpr (is_tree_)
        return wavelet_structure_tree();
      else
        return wavelet_structure_matrix();
    }

    auto ctx = ctx_t(size, levels);

    pc_in_external_parallel(text, size, levels, ctx);

    if constexpr (ctx_t::compute_zeros) {
      return wavelet_structure_matrix(std::move(ctx.bv()),
                                      std::move(ctx.zeros()));
    } else {
      return wavelet_structure_tree(std::move(ctx.bv()));
    }
  }
}; // class wx_pc

/******************************************************************************/
