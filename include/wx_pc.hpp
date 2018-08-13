/*******************************************************************************
 * include/wx_pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>

#include "util/ctx_single_level.hpp"
#include "util/pc.hpp"
#include "util/wavelet_structure.hpp"
#include "util/memory_types.hpp"

template <typename AlphabetType, bool is_tree_>
class wx_pc {

public:
  static constexpr bool  is_parallel = false;
  static constexpr bool  is_tree   = is_tree_;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);
  static constexpr bool  is_huffman_shaped = false;
  static constexpr memory_mode mem_mode = memory_mode::internal;

  template <typename InputType>
  static wavelet_structure compute(const InputType& text, const uint64_t size,
    const uint64_t levels) {

    using ctx_t = ctx_single_level<is_tree>;

    if(size == 0) {
      return wavelet_structure();
    }

    auto ctx = ctx_t(size, levels);

    pc(text, size, levels, ctx);

    if (ctx_t::compute_zeros) {
      return wavelet_structure(
        std::move(ctx.bv()), std::move(ctx.zeros()));
    } else {
      return wavelet_structure(std::move(ctx.bv()));
    }
  }
}; // class wx_pc

/******************************************************************************/
