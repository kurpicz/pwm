/*******************************************************************************
 * include/wx_ps_fe.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>

#include "arrays/memory_types.hpp"
#include "construction/ps_external.hpp"
#include "construction/wavelet_structure.hpp"

#include "wx_base.hpp"

template <typename AlphabetType, bool is_tree_, int word_packing_mode>
class wx_ps_fe : public wx_in_out_external<true, true> {
public:
  static constexpr bool is_parallel = false;
  static constexpr bool is_tree = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = false;

  template <typename InputType>
  static external_bit_vectors
  compute(const InputType& text, const uint64_t size, const uint64_t levels) {

    if (size == 0) {
      return external_bit_vectors();
    }
    return wx_ps_fe_builder<InputType, is_tree_, word_packing_mode>::build(
        text, size, levels);
  }
}; // class wx_ps

/******************************************************************************/
