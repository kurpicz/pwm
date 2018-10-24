/*******************************************************************************
 * include/queries/popcnt_traits.hpp
 *
 * Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 ******************************************************************************/

#pragma once

#include <limits>

struct popcnt_traits {
  static constexpr size_t bits_per_word = sizeof(uint64_t) * 8;
  static constexpr size_t l2_bit_size = 512;
  static constexpr size_t l1_bit_size = 4 * l2_bit_size;
  static constexpr size_t l0_bit_size = std::numeric_limits<uint32_t>::max();

  static constexpr size_t l2_block_cover = l2_bit_size / bits_per_word;
  static constexpr size_t l1_block_cover = l1_bit_size / bits_per_word;
  static constexpr size_t l0_block_cover = l0_bit_size / bits_per_word;
}; // struct popcnt_traits

/******************************************************************************/

