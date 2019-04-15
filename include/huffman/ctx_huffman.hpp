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
          typename borders_array,
          typename hist_array,
          template <bool> typename rho_type,
          bool require_initialization,
          template <bool> typename bv_type,
          typename huffman_codes>
class ctx_huffman {

public:
  ctx_huffman() = default;
  
  template<typename BvSizeArg>
  ctx_huffman(BvSizeArg&& size_arg,
              uint64_t const levels,
              rho_type<is_tree>&& rho,
              huffman_codes& codes,
              [[maybe_unused]] uint64_t const shards = 1)
    : hist_(levels + 1),
      borders_(levels + 1),
      zeros_(levels, 0),
      bv_(levels, size_arg),
      rho_(std::move(rho)),
      levels_(levels) {

    for (size_t i = 1; i < levels_; ++i) {
      for (auto const& cp : codes.code_pairs()) {
        if (cp.code_length >= i) {
          hist_.get(i)[cp.prefix(i)] = 0;
        }
      }
    }
  }

  uint64_t hist_size(uint64_t const level) {
    return 1ULL << level;
  }

  auto& hist_at_level(uint64_t const level) {
    return hist_.get(level);
  }

  auto hist_at_level(uint64_t const level) const {
    return hist_.get(level);
  }


  auto& borders_at_level(uint64_t const level) {
    return borders_.get(level);
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
    drop_me(std::move(borders_));
  }

  void discard_rho() {
    drop_me(std::move(rho_));
  }

private:
  hist_array hist_;
  borders_array borders_;
  std::vector<uint64_t>  zeros_;
  bv_type<require_initialization> bv_;
  rho_type<is_tree> rho_;
  uint64_t levels_;
}; // ctx_huffman

/******************************************************************************/
