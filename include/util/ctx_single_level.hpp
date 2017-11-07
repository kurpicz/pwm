/*******************************************************************************
 * include/util/ctx_single_level.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "bit_vectors.hpp"

// TODO: WM/WT abstract that selects zeros and rho
// TODO: flatten vectors where possible, to reduce indirection

// Overwrite information for each level
template<bool is_matrix>
struct ctx_single_level {
  std::vector<uint64_t> m_hist;
  std::vector<uint64_t> m_borders;
  std::vector<uint64_t> m_zeros;
  bit_vectors m_bv;
  std::vector<uint64_t> m_bit_reverse;

  ctx_single_level(uint64_t const size, uint64_t const levels)
  : m_hist(1ULL << levels, 0), m_borders(1ULL << levels, 0),
    m_zeros(levels, 0), m_bv(size, levels),
    m_bit_reverse(is_matrix ? BitReverse(levels - 1) : std::vector<uint64_t>(0)) { }

  uint64_t hist_size(uint64_t const level) {
    return 1ull << level;
  }

  uint64_t& hist(uint64_t const /*level*/, uint64_t const i) {
    return m_hist[i];
  }

  uint64_t hist(uint64_t const /*level*/, uint64_t const i) const {
    return m_hist[i];
  }

  uint64_t rho(size_t /*level*/, size_t i) {
    if (is_matrix) {
      return m_bit_reverse[i];
    }else {
      return i;
    }
  }

  void set_rho(size_t /*level*/, size_t i, uint64_t val) {
    m_bit_reverse[i] = val;
  }

  std::vector<uint64_t>& borders() {
    return m_borders;
  }

  static bool constexpr compute_zeros = is_matrix;
  static bool constexpr compute_rho = is_matrix;

  std::vector<uint64_t>& zeros() {
    return m_zeros;
  }

  bit_vectors& bv() {
    return m_bv;
  }

  bit_vectors const& bv() const {
    return m_bv;
  }

  void discard_non_merge_data() {
    // Not used in merge algorithm
  }
}; // ctx_single_level

/******************************************************************************/
