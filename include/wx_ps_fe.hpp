/*******************************************************************************
 * include/wx_ps_fe.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>
#include <sstream>
#include "construction/ps_external.hpp"
#include "construction/wavelet_structure_external.hpp"

#include "wx_base.hpp"

template <typename AlphabetType, bool is_tree_, int word_packing_mode>
class wx_ps_fe : public wx_in_out_external<true, true> {
public:
  static constexpr bool is_parallel = false;
  static constexpr bool is_tree = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = false;

  template <typename InputType>
  static wavelet_structure_external
  compute(const InputType& text,
          const uint64_t size,
          const uint64_t levels) {

    std::ostringstream name;
    name << "w" << (is_tree_ ? "t" : "m") << "_ps_fe_wp" << word_packing_mode;

    auto result = wavelet_structure_external::getNoHuffman(
        size, levels, is_tree_, 0, name.str());

    if (size > 0) {
      auto& bvs = wavelet_structure_external_writer::bvs(result);
      bvs.resize(0);
      wx_ps_fe_builder<InputType, is_tree_, word_packing_mode>::build(
          text, result);
    }

    return result;
  }
}; // class wx_ps

/******************************************************************************/
