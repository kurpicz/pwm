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
#include "permutation.hpp"

// TODO: WM/WT abstract that selects zeros and rho
// TODO: flatten vectors where possible, to reduce indirection

// Overwrite information for each level
template<bool is_matrix>
class ctx_single_level {

public:
  ctx_single_level(uint64_t const size, uint64_t const levels)
  : hist_(1ULL << levels, 0), borders_(1ULL << levels, 0),
    zeros_(levels, 0), bv_(size, levels),
    bit_reverse_(is_matrix ?
      bit_reverse_permutation(levels - 1) : std::vector<uint64_t>(0)) { }

  uint64_t hist_size(uint64_t const level) {
    return 1ull << level;
  }

  uint64_t& hist(uint64_t const /*level*/, uint64_t const i) {
    return hist_[i];
  }

  uint64_t hist(uint64_t const /*level*/, uint64_t const i) const {
    return hist_[i];
  }

  uint64_t rho(size_t /*level*/, size_t i) {
    if (is_matrix) {
      return bit_reverse_[i];
    }else {
      return i;
    }
  }

  void set_rho(size_t /*level*/, size_t i, uint64_t val) {
    bit_reverse_[i] = val;
  }

  std::vector<uint64_t>& borders() {
    return borders_;
  }

  static bool constexpr compute_zeros = is_matrix;
  static bool constexpr compute_rho = is_matrix;

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
    // Not used in merge algorithm
  }

private:
  std::vector<uint64_t> hist_;
  std::vector<uint64_t> borders_;
  std::vector<uint64_t> zeros_;
  bit_vectors bv_;
  std::vector<uint64_t> bit_reverse_;
}; // ctx_single_level

/******************************************************************************/
