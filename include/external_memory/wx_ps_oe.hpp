/*******************************************************************************
 * include/wx_ps_oe.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>
#include "external_memory/semi_out/ps_semi_external.hpp"
#include "external_memory/semi_out/inplace_ps_semi_external.hpp"
#include "wavelet_structure_external.hpp"

#include "wx_base.hpp"

template <typename AlphabetType, bool is_tree_, bool inplace_sort = false>
class wx_ps_oe : public wx_in_out_external<false, true, false, inplace_sort> {
public:
  static constexpr bool is_parallel = false;
  static constexpr bool is_tree = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = false;

  template <typename InputType>
  static wavelet_structure_external
  compute(const InputType& text, const uint64_t size, const uint64_t levels) {

    std::ostringstream name;
    name << "w" << (is_tree_ ? "t" : "m") << "_ps_oe";

    auto result =
        wavelet_structure_external_factory(is_tree_).
        histograms(true).
        construct(size, levels, name.str(), 0);

    if (size > 0) {
      if constexpr (inplace_sort)
        ps_out_external_inplace<AlphabetType, is_tree_>(text, result);
      else
        ps_out_external<AlphabetType, is_tree_>(text, result);
    }

    return result;
  }
}; // class wx_ps

/******************************************************************************/