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

template <typename AlphabetType, bool is_tree_, int word_packing_mode>
class wx_ps_fe {
public:

  WX_BASE(AlphabetType, is_tree_, false, false, memory_mode::external)

  template <typename InputType>
  static external_bit_vectors compute(const InputType& text, const uint64_t size,
    const uint64_t levels) {

    if(size == 0) { return external_bit_vectors(); }
    return wx_ps_fe_builder<InputType, is_tree_, word_packing_mode>::build(text, size, levels);
  }
}; // class wx_ps

/******************************************************************************/
