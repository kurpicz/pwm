/*******************************************************************************
 * include/wx_ps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>
#include <cmath>
#include "wx_base.hpp"
#include "util/ctx_single_level.hpp"
#include "util/wavelet_structure.hpp"
#include "util/ps.hpp"
#include "util/memory_types.hpp"


template <typename AlphabetType, bool is_tree_>
class wx_ps {
public:

  WX_BASE(AlphabetType, is_tree_, false, false, memory_mode::internal)

  template <typename InputType>
  static wavelet_structure compute(const InputType& text, const uint64_t size,
    const uint64_t levels) {

    using ctx_t = ctx_single_level<is_tree_>;

    if(size == 0) { return wavelet_structure(); }

    auto ctx = ctx_t(size, levels);

    auto sorted_text = std::vector<AlphabetType>(size);
    ps(text, size, levels, ctx, sorted_text.data());

    if (ctx_t::compute_zeros)  {
      return wavelet_structure(std::move(ctx.bv()), std::move(ctx.zeros()));
    } else {
      return wavelet_structure(std::move(ctx.bv()));
    }
  }
}; // class wx_ps

/******************************************************************************/
