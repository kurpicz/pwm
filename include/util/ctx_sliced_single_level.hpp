/*******************************************************************************
 * include/util/ctx_sliced_single_level.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "common.hpp"

// TODO: WM/WT abstract that selects zeros and rho
// TODO: flatten vectors where possible, to reduce indirection

template <bool is_matrix>
class ctx_sliced_single_level {
public:
  ctx_sliced_single_level() = default;

  ctx_sliced_single_level(uint64_t const size, uint64_t const levels,
    const uint64_t omp_size)
    : hist_(omp_size, std::vector<uint64_t>(1ULL << levels, 0)),
      borders_(omp_size, std::vector<uint64_t>(1ULL << levels, 0)),
      zeros_(levels, 0), bv_(size, levels),
      bit_reverse_(
      is_matrix ? BitReverse(levels - 1) : std::vector<uint64_t>(0)) { }

  uint64_t borders_size(uint64_t const level) {
    return 1ULL << level;
  }

  uint64_t& borders(uint64_t const rank, uint64_t const index) {
    return borders_[rank][index];
  }

  uint64_t borders(uint64_t const rank, uint64_t const index) const {
    return borders_[rank][index];
  }

  uint64_t hist_size(uint64_t const level) {
    return 1ULL << level;
  }

  uint64_t& hist(uint64_t const rank, uint64_t const index) {
    return hist_[rank][index];
  }

  uint64_t hist(uint64_t const rank, uint64_t const index) const {
    return hist_[rank][index];
  }

  uint64_t rho (uint64_t /*level*/, uint64_t const index) {
    if (is_matrix) {
      return bit_reverse_[index];
    } else {
      return index;
    }
  }

  void set_rho(uint64_t /*level*/, uint64_t const index, uint64_t const val) {
    bit_reverse_[index] = val;
  }

  static bool constexpr compute_zeros = is_matrix;
  static bool constexpr compute_rho = is_matrix;

  std::vector<uint64_t>& zeros() {
    return zeros_;
  }

  Bvs& bv() {
    return bv_;
  }

  Bvs const& bv() const {
    return bv_;
  }

  void discard_non_merge_data() {
    // Not used in merge algorithms.
  }

private:
  std::vector<std::vector<uint64_t>> hist_;
  std::vector<std::vector<uint64_t>> borders_;

  std::vector<uint64_t> zeros_;
  Bvs bv_;
  std::vector<uint64_t> bit_reverse_;
}; // class ctx_sliced_single_level 

/******************************************************************************/
