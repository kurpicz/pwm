/*******************************************************************************
 * include/util/ctx_compute_borders.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "bit_vectors.hpp"

// TODO: WM/WT abstract that selects zeros and rho
// TODO: flatten vectors where possible, to reduce indirection

/// Keep calculated information for individual levels around
template<bool is_matrix>
class ctx_compute_borders {

public:
  using rho_t = typename rho_dispatch<is_matrix>::type;

  ctx_compute_borders() = default;

  ctx_compute_borders(uint64_t const size, uint64_t const levels,
    rho_t const& rho)
  : hist_(levels + 1, std::vector<uint64_t>(2)), rho_(&rho),
    borders_(levels + 1, std::vector<uint64_t>(2)),
    zeros_(levels, 0), bv_(levels, size), levels_(levels) {

    for(uint64_t level = 0; level < levels_ + 1; ++level) {
      hist_[level].reserve(hist_size(level));
      hist_[level].resize(hist_size(level));
      borders_[level].reserve(hist_size(level));
      borders_[level].resize(hist_size(level));
    }
  }

  static bool constexpr compute_zeros = is_matrix;
  static bool constexpr compute_rho = false;

  void fill_borders() {
    for (uint64_t level = levels_ - 1; level > 0; --level) {
      for (uint64_t pos = 0; pos < hist_size(level); ++pos) {
        hist_[level][pos] = hist_[level + 1][pos << 1] +
          hist_[level + 1][(pos << 1) + 1];
      }

      borders_[level][0] = 0;
      for (uint64_t pos = 1; pos < hist_size(level); ++pos) {
        auto const prev_rho = rho(level, pos - 1);

        borders_[level][rho(level, pos)] =
          borders_[level][prev_rho] + hist_[level][prev_rho];
      }

      // The number of 0s is the position of the first 1 in the previous level
      if (compute_zeros) {
        zeros_[level - 1] = borders_[level][1];
      }
    }
  }

  inline uint64_t hist_size(uint64_t const level) {
    return 1ULL << level;
  }

  uint64_t& hist(uint64_t const level, uint64_t const i) {
    return hist_[level][i];
  }

  uint64_t hist(uint64_t const level, uint64_t const i) const {
    return hist_[level][i];
  }

  uint64_t rho(uint64_t level, uint64_t i) {
    return (*rho_)(level, i);
  }

  void set_rho(uint64_t /*level*/, uint64_t /*i*/, uint64_t /*val*/) {
    // rho_ is already calculated
  }

  uint64_t& borders(uint64_t const level, uint64_t const i) {
    return borders_[level][i];
  }

  uint64_t borders(uint64_t const level, uint64_t const i) const {
    return borders_[level][i];
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
  std::vector<std::vector<uint64_t>> hist_;
  rho_t const* rho_ = nullptr;
  std::vector<std::vector<uint64_t>> borders_;
  std::vector<uint64_t> zeros_;
  bit_vectors bv_;
  uint64_t levels_;

}; // ctx_compute_borders

/******************************************************************************/
