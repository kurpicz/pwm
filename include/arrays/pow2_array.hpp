/*******************************************************************************
 * include/arrays/pow2_array.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "flat_two_dim_array.hpp"
#include "util/common.hpp"

struct pow2_level_size {
  static uint64_t level_size(const uint64_t level) {
    return 1ULL << level;
  };

  static constexpr bool is_bit_vector = false;
}; // struct pow2_level_size

using pow2_array = flat_two_dim_array<uint64_t, pow2_level_size>;

/******************************************************************************/
