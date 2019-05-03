/*******************************************************************************
 * include/util/ctx_all_levels.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "arrays/bit_vectors.hpp"
#include "arrays/pow2_array.hpp"
#include "util/permutation.hpp"

// TODO: WM/WT abstract that selects zeros and rho

// Context for partial results:
//    Keep histograms of all levels
//    Keep borders of current level only
//    Make borders word aligned
//    Only tree implementation
class ctx_partial {
private:
  uint64_t size_;
  const uint64_t levels_;
  const uint64_t max_overhead_; // 63 bits per WT interval

  pow2_array hist_;
  std::vector<uint64_t> borders_;
  bit_vectors bv_;

public:

  ctx_partial() = default;

  ctx_partial(uint64_t const size, uint64_t const levels)
      : size_(size), levels_(levels),
        max_overhead_(((1ULL << levels) - 1) * 63),
        hist_(levels_ + 1),
        borders_(1ULL << levels_, 0),
        bv_(levels_, size_ + max_overhead_) {}

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

  pow2_array extract_final_hist() {
    return std::move(hist_);
  }

  uint64_t size() const {
    return size_;
  }

  uint64_t levels() const {
    return levels_;
  }
}; // class ctx_all_levels

/******************************************************************************/
