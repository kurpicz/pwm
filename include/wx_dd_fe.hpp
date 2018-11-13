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

template <typename AlphabetType, bool is_tree_>
class wx_dd_fe : public wx_in_out_external<true, true> {
public:
  static constexpr bool is_parallel = true;
  static constexpr bool is_tree = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = false;

  template <typename InputType>
  static wavelet_structure_external
  compute(const InputType& text,
          const uint64_t size,
          const uint64_t levels) {

    std::ostringstream name;
    name << "w" << (is_tree_ ? "t" : "m") << "_dd_fe";

    auto result =
        wavelet_structure_external_factory(is_tree_).
            histograms().zeros().
            construct(size, levels, name.str(), 0);

    if (size > 0) {
      //TODO: implement
    }

    return result;
  }
}; // class wx_ps

/******************************************************************************/
