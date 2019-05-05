/*******************************************************************************
 * include/construction/pc_dd_fe/ctx_partial.hpp
 *
 * Copyright (C) 2019 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "arrays/bit_vectors.hpp"
#include "arrays/pow2_array.hpp"

// Context for partial results:
//    Keep histograms of all levels
//    Keep borders of current level only
//    Make borders word aligned
//    Only tree implementation
class ctx_partial {
private:
  const uint64_t size_;
  const uint64_t levels_;
  const uint64_t max_overhead_; // 63 bits per WT interval

  pow2_array hist_;
  std::vector<uint64_t> borders_;
  bit_vectors bv_;

  uint64_t data_size_;
  std::vector<uint64_t> level_data_sizes_;

public:

  ctx_partial() = default;

  ctx_partial(uint64_t const size, uint64_t const levels)
      : size_(size), levels_(levels),
        max_overhead_(((1ULL << levels) - 1) * 63),
        hist_(levels_ + 1),
        borders_(1ULL << levels_, 0),
        bv_(levels_, size_ + max_overhead_),
        data_size_(0),
        level_data_sizes_(levels, 0) {}

  uint64_t hist_size(uint64_t const level) {
    return 1ULL << level;
  }

  uint64_t& hist(uint64_t const level, uint64_t const i) {
    return hist_[level][i];
  }

  uint64_t hist(uint64_t const level, uint64_t const i) const {
    return hist_[level][i];
  }

  std::vector<uint64_t>& borders() {
    return borders_;
  }

  bit_vectors& bv() {
    return bv_;
  }

  bit_vectors const& bv() const {
    return bv_;
  }

  pow2_array& hist() {
    return hist_;
  }

  std::vector<uint64_t> flat_hist() {
    std::vector<uint64_t> result((1ULL << levels_) - 1);
    uint64_t i = 0;
    for (uint64_t l = 0; l < levels_; ++l) {
      const uint64_t level_hist_size = hist_size(l);
      for (uint64_t h = 0; h < level_hist_size; ++h) {
        result[i++] = hist_[l][h];
      }
    }
    return result;
  }

//  pow2_array extract_final_hist() {
//    return std::move(hist_);
//  }

  uint64_t& data_size() {
    return data_size_;
  }

  uint64_t data_size() const {
    return data_size_;
  }

  std::vector<uint64_t>& level_data_sizes() {
    return level_data_sizes_;
  }

  const std::vector<uint64_t>& level_data_sizes() const {
    return level_data_sizes_;
  }

  uint64_t size() const {
    return size_;
  }

  uint64_t levels() const {
    return levels_;
  }
}; // class ctx_all_levels

/******************************************************************************/
