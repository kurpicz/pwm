/*******************************************************************************
 * include/util/debug.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <iostream>
#include <string>

#include "util/wavelet_structure.hpp"
#include "util/debug.hpp"

[[gnu::unused]] // TODO: C++17 [[maybe_unused]]
static void print_structure(std::ostream& out, wavelet_structure const& structure) {
  const base_bit_vectors& bv = structure.bvs();
  const std::vector<uint64_t>& zeros = structure.zeros();
  for (uint64_t i = 0; i < bv.levels(); i++) {
    out << "   bv["<<i<<"]";

    out << "[";
    for (uint64_t j = 0; j < bv.level_bit_size(i); j++) {
      out << uint64_t(bit_at(bv[i], j)) << "";
    }
    out << "]";
    if (!structure.is_tree()) {
        out << " zeros[" << i << "] = " << zeros.at(i);
    }
    out << std::endl;
  }
  out << std::endl;
}

/******************************************************************************/
