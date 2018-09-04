/*******************************************************************************
 * include/wx_ps_fe.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>
#include "wx_base.hpp"
#include "util/ctx_single_level_external.hpp"
#include "util/wavelet_structure.hpp"
#include "util/ps_external.hpp"
#include "util/memory_types.hpp"

template <typename AlphabetType, bool is_tree_>
class wx_ps_fe {
public:

  WX_BASE(AlphabetType, is_tree_, false, false, memory_mode::external)

  template <typename InputType>
  static external_bit_vectors compute(const InputType& text, const uint64_t size,
    const uint64_t levels) {

    using ctx_t = ctx_single_level_external<is_tree_>;

    if(size == 0) { return external_bit_vectors(); }

    auto ctx = ctx_t(size, levels);
    return ps_fully_external2<AlphabetType>(text, size, levels, ctx);
  }
}; // class wx_ps

/******************************************************************************/
