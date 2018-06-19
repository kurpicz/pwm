/*******************************************************************************
 * include/util/common.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <stdint.h>
#include <vector>
#include <cassert>
#include <type_traits>
#include <climits>

#include "permutation.hpp"

constexpr uint64_t word_size(uint64_t size) {
  return (size + 63ULL) >> 6;
}

// A little template helper for dropping a type early
template <typename T>
void drop_me(T const&) = delete;
template <typename T>
void drop_me(T&) = delete;
template <typename T>
void drop_me(T const&&) = delete;
template <typename T>
void drop_me(T&& t) {
  std::remove_reference_t<T>(std::move(t));
}

// This is not compatible with STXXL and has been replaced with pwm_log2
//~ constexpr uint64_t log2(uint64_t n) {
  //~ return (n < 2) ? 1 : 1 + log2(n / 2);
//~ }

constexpr uint64_t pwm_log2(uint64_t n) {
  return (n < 2) ? 1 : 1 + pwm_log2(n / 2);
}

template<typename WordType = uint64_t, typename bv_t>
inline auto bit_at(const bv_t& bv, uint64_t i) -> bool {
  constexpr WordType BITS = (sizeof(WordType) * CHAR_BIT);
  constexpr WordType MOD_MASK = BITS - 1;

  uint64_t offset = i / BITS;
  uint64_t word_offset = i & MOD_MASK;
  return (bv[offset] >> (MOD_MASK - word_offset)) & 1ULL;
}

/******************************************************************************/
