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

template <bool is_tree,
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
    //TODO Init umaps with zeros at required places

    for (auto const& cp : codes.code_pairs()) {
      for (uint64_t i = 1; i < std::min(cp.code_length, levels_); ++i) {
        hist_[i][cp.prefix(i)] = 0;
        borders_[i][cp.prefix(i)] = 0;
      }
    }
  }

  uint64_t hist_size(uint64_t const level) {
    return 1ULL << level;
  }

  auto& hist_at_level(uint64_t const level) {
    return hist_[level];
  }

  auto hist_at_level(uint64_t const level) const {
    return hist_[level];
  }


  auto& borders_at_level(uint64_t const level) {
    return borders_[level];
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
  using help_array = std::vector<std::unordered_map<uint64_t, uint64_t>>;

  help_array hist_;
  help_array borders_;
  std::vector<uint64_t>  zeros_;
  bv_type<require_initialization> bv_;
  rho_type<is_tree> rho_;
  uint64_t levels_;
}; // ctx_huffman

/******************************************************************************/
