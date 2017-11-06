/*******************************************************************************
 * include/util/ctx_all_levels.hpp
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

/// Keep calculated information for individual levels around
template<bool is_matrix, typename rho_t>
struct ctx_all_levels {
  std::vector<std::vector<uint64_t>> m_hist;
  rho_t const* m_rho = nullptr;
  std::vector<uint64_t> m_borders;
  std::vector<uint64_t> m_zeros;
  bit_vectors m_bv;

  ctx_all_levels() = default;

  ctx_all_levels(uint64_t const size, uint64_t const levels, rho_t const& rho)
  : m_hist(levels + 1, std::vector<uint64_t>(2)),
    m_rho(&rho), m_borders(1ULL << levels, 0), m_zeros(levels, 0),
    m_bv(size, levels) {

    for(size_t level = 0; level < (levels + 1); level++) {
      m_hist[level].reserve(hist_size(level));
      m_hist[level].resize(hist_size(level));
    }
  }

  uint64_t hist_size(uint64_t const level) {
    return 1ull << level;
  }

  uint64_t& hist(uint64_t const level, uint64_t const i) {
    return m_hist[level][i];
  }

  uint64_t hist(uint64_t const level, uint64_t const i) const {
    return m_hist[level][i];
  }

  uint64_t rho(size_t level, size_t i) {
    return (*m_rho)(level, i);
  }

  void set_rho(size_t /*level*/, size_t /*i*/, uint64_t /*val*/) {
    // m_rho is already calculated
  }

  std::vector<uint64_t>& borders() {
    return m_borders;
  }

  static bool constexpr compute_zeros = is_matrix;
  static bool constexpr compute_rho = false;

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
    drop_me(std::move(m_borders));
  }
}; // ctx_all_levels

/******************************************************************************/
