/*******************************************************************************
 * include/wx_pc_ss.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/wavelet_structure.hpp"

template <typename AlphabteType, bool is_matrix>
class wx_pc_ss {

public:
  static constexpr bool is_parallel = false;
  static constexpr bool is_tree     = !is_matrix;
  static constexpr uint8_t word_width = sizeof(AlphabteType);

  static wavelet_structure compute(AlphabteType const* const text,
                                   const uint64_t size,
                                   const uint64_t levels) {
    if (size == 0) {
      return wavelet_structure();
    }

    auto ctx = ctx_t(size, levels);
  }
}; // class wc_pc_ss
