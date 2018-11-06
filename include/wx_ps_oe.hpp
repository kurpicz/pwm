/*******************************************************************************
 * include/wx_ps_oe.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <bitset>
#include <vector>

#include "arrays/memory_types.hpp"
#include "construction/ctx_single_level.hpp"
#include "construction/ps_external.hpp"
#include "construction/wavelet_structure.hpp"

#include "wx_base.hpp"

template <typename AlphabetType, bool is_tree_>
class wx_ps_oe : public wx_in_out_external<false, true> {
public:
  static constexpr bool is_parallel = false;
  static constexpr bool is_tree = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = false;

  template <typename InputType>
  static external_bit_vectors
  compute(const InputType& text, const uint64_t size, const uint64_t levels) {

    using ctx_t = ctx_single_level<is_tree_>;

    if (size == 0) {
      return external_bit_vectors();
    }

    auto ctx = ctx_t(size, levels);
    return ps_out_external<AlphabetType>(text, size, levels, ctx);
  }
}; // class wx_ps

/******************************************************************************/
