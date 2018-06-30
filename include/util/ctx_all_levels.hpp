/*******************************************************************************
 * include/util/ctx_all_levels.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "permutation.hpp"
#include "bit_vectors.hpp"
#include "border_hist_array.hpp"

// TODO: WM/WT abstract that selects zeros and rho

/// Keep calculated information for individual levels around
template<typename OutputType, bool is_tree>
class ctx_all_levels {
  
public:
  using bit_vectors = OutputType;
  using rho_t = typename rho_dispatch<is_tree>::type;

  ctx_all_levels() = default;

  ctx_all_levels(uint64_t const size, uint64_t const levels, rho_t const& rho)
  : hist_(levels + 1), rho_(&rho),
    borders_(1ULL << levels, 0), zeros_(levels, 0), bv_(levels, size) { }

  static bool constexpr compute_zeros = !is_tree;
  static bool constexpr compute_rho = false;

  uint64_t hist_size(uint64_t const level) {
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

  bit_vectors& bv() {
    return bv_;
  }

  bit_vectors const& bv() const {
    return bv_;
  }

  void discard_non_merge_data() {
    drop_me(std::move(borders_));
  }

private:
  border_hist_array hist_;
  rho_t const* rho_ = nullptr;
  std::vector<uint64_t> borders_;
  std::vector<uint64_t> zeros_;
  bit_vectors bv_;
}; // class ctx_all_levels

/******************************************************************************/
