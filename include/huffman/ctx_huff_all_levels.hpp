/*******************************************************************************
 * include/huffman/ctx_huff_all_levels.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "arrays/bit_vectors.hpp"
#include "arrays/pow2_array.hpp"
#include "huff_bit_vectors.hpp"

/// Keep calculated information for individual levels around
template <bool is_tree>
class ctx_huff_all_levels {

public:
  using rho_t = typename rho_dispatch<is_tree>::type;

  ctx_huff_all_levels() = default;

  ctx_huff_all_levels(const std::vector<uint64_t>& sizes,
                      uint64_t const levels,
                      rho_t const& rho)
      : hist_(levels + 1),
        rho_(&rho),
        borders_(1ULL << levels, 0),
        zeros_(levels, 0),
        bv_(levels, sizes) {}

  static bool constexpr compute_zeros = !is_tree;
  static bool constexpr compute_rho = false;

  uint64_t hist_size(uint64_t const level) const {
    return 1ULL << level;
  }

  uint64_t& hist(uint64_t const level, uint64_t const i) {
    return hist_[level][i];
  }

  uint64_t hist(uint64_t const level, uint64_t const i) const {
    return hist_[level][i];
  }

  uint64_t rho(size_t level, size_t i) {
    return (*rho_)(level, i);
  }

  void set_rho(size_t /*level*/, size_t /*i*/, uint64_t /*val*/) {
    // rho_ is already calculated
  }

  std::vector<uint64_t>& borders() {
    return borders_;
  }

  std::vector<uint64_t>& zeros() {
    return zeros_;
  }

  huff_bit_vectors& bv() {
    return bv_;
  }

  huff_bit_vectors const& bv() const {
    return bv_;
  }

  void discard_non_merge_data() {
    drop_me(std::move(borders_));
  }

private:
  pow2_array hist_;
  rho_t const* rho_ = nullptr;
  std::vector<uint64_t> borders_;
  std::vector<uint64_t> zeros_;
  huff_bit_vectors bv_;
}; // class ctx_huff_all_levels

/******************************************************************************/
