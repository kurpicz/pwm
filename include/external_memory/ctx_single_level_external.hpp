/*******************************************************************************
 * include/external_memory/ctx_single_level_external.hpp
 *
 * Copyright (C) 2019 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/permutation.hpp"

template <bool is_tree>
class ctx_single_level_external {
private:
  std::vector<uint64_t> borders_;
  std::vector<uint64_t> bit_reverse_;
public:
  ctx_single_level_external(uint64_t /*size*/, uint64_t levels)
      : borders_(1ULL << levels, 0),
        bit_reverse_(is_tree ? std::vector<uint64_t>(0)
                             : bit_reverse_permutation(levels - 1)) {

  }

  uint64_t rho(size_t /*level*/, size_t i) const {
    if constexpr (!is_tree) {
      return bit_reverse_[i];
    } else {
      return i;
    }
  }

  void set_rho(size_t /*level*/, size_t i, uint64_t val) {
    if constexpr (!is_tree) {
      bit_reverse_[i] = val;
    }
  }

  std::vector<uint64_t>& borders() {
    return borders_;
  }
};

/******************************************************************************/
