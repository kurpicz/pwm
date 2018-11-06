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
#include "construction/wavelet_structure.hpp"

#include "wx_base.hpp"

template <typename AlphabetType, typename ParallelAlgorithmType>
class wx_dd_fe : public wx_in_out_external<true, true> {
public:
  static constexpr bool is_parallel = ParallelAlgorithmType::is_parallel;
  static constexpr bool is_tree = ParallelAlgorithmType::is_tree;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = false;

  template <typename InputType>
  static external_bit_vectors
  compute(const InputType& text, const uint64_t size, const uint64_t levels) {

    if (size == 0) {
      return external_bit_vectors();
    }
    //TODO:
    return external_bit_vectors();
  }
}; // class wx_ps

/******************************************************************************/
