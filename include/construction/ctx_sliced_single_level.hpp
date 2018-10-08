/*******************************************************************************
 * include/util/ctx_sliced_single_level.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "arrays/bit_vectors.hpp"
#include "arrays/helper_array.hpp"
#include "util/permutation.hpp"

// TODO: WM/WT abstract that selects zeros and rho

template <bool is_tree>
class ctx_sliced_single_level {

public:
  ctx_sliced_single_level() = default;

  ctx_sliced_single_level(uint64_t const size,
                          uint64_t const levels,
                          const uint64_t omp_size)
      : hist_(omp_size, 1ULL << levels),
        borders_(omp_size, 1ULL << levels),
        zeros_(levels, 0),
        bv_(levels, size),
        bit_reverse_(is_tree ? std::vector<uint64_t>(0)
                             : bit_reverse_permutation(levels - 1)) {}

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

  uint64_t rho(uint64_t /*level*/, uint64_t const index) {
    if constexpr (!is_tree) {
      return bit_reverse_[index];
    } else {
      return index;
    }
  }

  void set_rho(uint64_t /*level*/, uint64_t const index, uint64_t const val) {
    bit_reverse_[index] = val;
  }

  static bool constexpr compute_zeros = !is_tree;
  static bool constexpr compute_rho = !is_tree;

  std::vector<uint64_t>& zeros() {
    return zeros_;
  }

  bit_vectors& bv() {
    return bv_;
  }

  bit_vectors const& bv() const {
    return bv_;
  }

  void discard_non_merge_data() {
    // Not used in merge algorithms.
  }

private:
  helper_array hist_;
  helper_array borders_;
  std::vector<uint64_t> zeros_;
  bit_vectors bv_;
  std::vector<uint64_t> bit_reverse_;
}; // class ctx_sliced_single_level

/******************************************************************************/
