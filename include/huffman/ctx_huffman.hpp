/*******************************************************************************
 * include/huffman/ctx_huffman.hpp
 *
 * Copyright (C) 2019 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <unordered_map>
#include <vector>

#include "util/debug_assert.hpp"

namespace ctx_huffman_options{

struct single_level : public std::unordered_map<uint64_t, uint64_t> {
  inline single_level() = default;

  inline single_level(uint64_t /*levels*/) { }

  inline auto& get(uint64_t /*level*/) {
    return *this;
  }

  inline auto& get(uint64_t /*level*/) const {
    return *this;
  }

  static constexpr bool is_not_sharded = true;
  static constexpr bool is_not_leveled = true;
};

struct all_level : public std::vector<std::unordered_map<uint64_t, uint64_t>> {
  inline all_level() = default;

  inline all_level(uint64_t const levels)
    : std::vector<std::unordered_map<uint64_t, uint64_t>>(levels) { }

  inline auto& get(uint64_t level) {
    return (*this)[level];
  }
  inline auto& get(uint64_t level) const {
    return (*this)[level];
  }
  static constexpr bool is_not_sharded = true;
  static constexpr bool is_not_leveled = false;
};

namespace huff_borders {
  using ctx_huffman_options::single_level;
  using ctx_huffman_options::all_level;
}

namespace huff_hist {
  using ctx_huffman_options::all_level;
}

} // namespace ctx_huffman_options

template <bool is_tree,
          typename upper_borders_array,
          typename lower_borders_array,
          typename hist_array,
          template <bool> typename rho_type,
          bool require_initialization,
          template <bool> typename bv_type,
          typename huffman_codes>
class ctx_huffman {

public:
  // Number of levels that we use consecutive memory to store the histogram and
  // borders, because access is faster.
  static constexpr uint64_t MAX_UPPER_LEVELS = 7;
  static constexpr uint64_t UPPER_LEVEL_SIZE = MAX_UPPER_LEVELS + 1;

  ctx_huffman() = default;
  
  template<typename BvSizeArg>
  ctx_huffman(BvSizeArg&& size_arg,
              uint64_t const levels,
              rho_type<is_tree>&& rho,
              huffman_codes& codes,
              [[maybe_unused]] uint64_t const shards = 1)
    : hist_(levels + 1),
      upper_borders_(std::min(levels, UPPER_LEVEL_SIZE) + 1, 1),
      lower_borders_((levels > UPPER_LEVEL_SIZE ?
                      levels - UPPER_LEVEL_SIZE : 0) + 1),
      zeros_(levels, 0),
      bv_(levels, size_arg),
      rho_(std::move(rho)),
      levels_(levels) {

    for (size_t i = 1; i < levels_; ++i) {
      for (auto const& cp : codes.code_pairs()) {
        if (cp.code_length() >= i) {
          hist_.get(i)[cp.prefix(i)] = 0;
        }
      }
    }
  }

  uint64_t hist_size(uint64_t const level) const {
    return 1ULL << level;
  }

  auto& hist_at_level(uint64_t const level) {
    return hist_.get(level);
  }

  auto hist_at_level(uint64_t const level) const {
    return hist_.get(level);
  }

  auto& lower_borders_at_level(uint64_t const level) {
    DCHECK_GT(level, MAX_UPPER_LEVELS);
    return lower_borders_.get(level - MAX_UPPER_LEVELS);
  }

  auto upper_borders_at_level(uint64_t const level) {
    DCHECK_LE(level, MAX_UPPER_LEVELS);
    return upper_borders_.get(0 /* TODO: 0 is arbitrary here*/, level);
  }


  uint64_t rho(size_t level, size_t i) {
    return rho_.get(level, i);
  }

  void set_rho(size_t level, size_t i, uint64_t val) {
    rho_.set(level, i, val);
  }

  static bool constexpr compute_zeros = !is_tree;
  static bool constexpr compute_rho = rho_type<is_tree>::compute_rho;

  auto zeros() {
    return span<uint64_t>(zeros_);
  }

  auto zeros() const {
    return span<uint64_t const>(zeros_);
  }

  std::vector<uint64_t>&& take_zeros() {
    return std::move(zeros_);
  }

  bv_type<require_initialization>& bv() {
    return bv_;
  }

  bv_type<require_initialization> const& bv() const {
    return bv_;
  }

  void discard_borders() {
    drop_me(std::move(upper_borders_));
    drop_me(std::move(lower_borders_));
  }

  void discard_rho() {
    drop_me(std::move(rho_));
  }

private:
  hist_array hist_;
  upper_borders_array upper_borders_;
  lower_borders_array lower_borders_;
  std::vector<uint64_t>  zeros_;
  bv_type<require_initialization> bv_;
  rho_type<is_tree> rho_;
  uint64_t levels_;
}; // ctx_huffman

/******************************************************************************/
