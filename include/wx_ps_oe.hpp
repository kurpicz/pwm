/*******************************************************************************
 * include/wx_ps_oe.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>
#include <bitset>
#include "wx_base.hpp"
#include "util/ctx_single_level.hpp"
#include "util/wavelet_structure.hpp"
#include "util/ps_external.hpp"
#include "util/memory_types.hpp"

template <typename AlphabetType, bool is_tree_>
class wx_ps_oe {
public:

  WX_BASE(AlphabetType, is_tree_, false, false, memory_mode::external_output)

  template <typename InputType>
  static external_bit_vectors compute(const InputType& text, const uint64_t size,
    const uint64_t levels) {

    using ctx_t = ctx_single_level<is_tree_>;

    if(size == 0) { return external_bit_vectors(); }

    auto ctx = ctx_t(size, levels);

    auto sorted_text = std::vector<AlphabetType>(size);
    return ps_out_external(text, size, levels, ctx, sorted_text.data());
  }
}; // class wx_ps

/******************************************************************************/
