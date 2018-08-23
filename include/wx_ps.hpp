/*******************************************************************************
 * include/wx_ps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>

#include "util/ctx_single_level.hpp"
#include "util/wavelet_structure.hpp"
#include "util/ps.hpp"

template <typename AlphabetType, bool is_tree_>
class wx_ps {

public:
  static constexpr bool  is_parallel = false;
  static constexpr bool  is_tree   = is_tree_;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);
  static constexpr bool  is_huffman_shaped = false;

  using ctx_t = ctx_single_level<is_tree>;

  template <typename InputType>
  static wavelet_structure compute(const InputType& text, const uint64_t size,
    const uint64_t levels) {

    if(size == 0) { return wavelet_structure(); }

    auto ctx = ctx_t(size, levels);

    auto sorted_text = std::vector<AlphabetType>(size);
    ps(text, size, levels, ctx, sorted_text.data());

    if (ctx_t::compute_zeros)  {
      return wavelet_structure(std::move(ctx.bv()), std::move(ctx.zeros()), false);
    } else {
      return wavelet_structure(std::move(ctx.bv()), false);
    }
  }
}; // class wx_ps

/******************************************************************************/
