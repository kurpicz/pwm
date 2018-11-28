/*******************************************************************************
 * include/util/ctx_generic.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <type_traits>

#include "arrays/bit_vectors.hpp"
#include "arrays/helper_array.hpp"
#include "arrays/pow2_array.hpp"

#include "util/permutation.hpp"

namespace ctx_options {
template <bool is_tree>
class pre_computed_rho {
  using rho_t = typename rho_dispatch<is_tree>::type;

  rho_t const* rho_ = nullptr;
public:
  inline uint64_t get(size_t level, size_t i) const {
    return (*rho_)(level, i);
  }

  inline void set(size_t /*level*/, size_t /*i*/, uint64_t /*value*/) {
    DCHECK(false); // Should not be called
  }

  inline pre_computed_rho() = default;
  inline pre_computed_rho(rho_t const& rho) : rho_(&rho) {}

  static bool constexpr compute_rho = false;
};

template <bool is_tree>
class live_computed_rho {
  std::vector<uint64_t> bit_reverse_;
public:
  inline uint64_t get(size_t /*level*/, size_t i) const {
    if constexpr (!is_tree) {
      return bit_reverse_[i];
    } else {
      return i;
    }
  }

  inline void set(size_t /*level*/, [[maybe_unused]] size_t i, [[maybe_unused]] uint64_t value) {
    if constexpr (!is_tree) {
      bit_reverse_[i] = value;
    } else {
      DCHECK(false); // Should not be called
    }
  }

  inline live_computed_rho() = default;
  inline live_computed_rho(uint64_t levels)
      : bit_reverse_(is_tree ? std::vector<uint64_t>(0)
                             : bit_reverse_permutation(levels - 1)) {}

  static bool constexpr compute_rho = !is_tree;
};

struct single_level : public std::vector<uint64_t> {
  inline single_level() = default;
  inline single_level(uint64_t const levels, uint64_t const /*shards*/)
      : std::vector<uint64_t>(1ULL << levels, 0) {}
  inline auto get(uint64_t /*shard*/, uint64_t /*level*/) {
    return span<uint64_t>(this->data(), this->size());
  }
  inline auto get(uint64_t /*shard*/, uint64_t /*level*/) const {
    return span<uint64_t const>(this->data(), this->size());
  }
  static constexpr bool is_not_sharded = true;
  static constexpr bool is_not_leveled = true;
};
struct all_level : public pow2_array {
  inline all_level() = default;

  inline all_level(uint64_t const levels, uint64_t const /*shards*/)
      : pow2_array(levels + 1) {}
  inline auto get(uint64_t /*shard*/, uint64_t level) {
    return (*this)[level];
  }
  inline auto get(uint64_t /*shard*/, uint64_t level) const {
    return (*this)[level];
  }
  static constexpr bool is_not_sharded = true;
  static constexpr bool is_not_leveled = false;
};
struct sharded_single_level : public helper_array {
  inline sharded_single_level() = default;

  inline sharded_single_level(uint64_t const levels, uint64_t const shards)
      : helper_array(shards, 1ULL << levels) {}
  inline auto get(uint64_t shard, uint64_t /*level*/) {
    return (*this)[shard];
  }
  inline auto get(uint64_t shard, uint64_t /*level*/) const {
    return (*this)[shard];
  }
  static constexpr bool is_not_sharded = false;
  static constexpr bool is_not_leveled = true;
};

namespace borders {
    using ctx_options::single_level;
    using ctx_options::all_level;
    using ctx_options::sharded_single_level;
}
namespace hist {
    using ctx_options::single_level;
    using ctx_options::all_level;
    using ctx_options::sharded_single_level;
}
constexpr bool bv_initialized = true;
constexpr bool bv_uninitialized = false;

} // namespace ctx_options

template <bool is_tree,
          typename borders_array,
          typename hist_array,
          template <bool> typename rho_type,
          bool require_initialization,
          template <bool> typename bv_type>
class ctx_generic {
public:
  ctx_generic() = default;

  template<typename BvSizeArg>
  ctx_generic(BvSizeArg&& size_arg,
              uint64_t const levels,
              rho_type<is_tree>&& rho,
              uint64_t const shards = 1)
      : hist_(levels, shards),
        borders_(levels, shards),
        zeros_(levels, 0),
        bv_(levels, size_arg),
        rho_(std::move(rho)),
        levels_(levels) {}

  uint64_t hist_size(uint64_t const level) {
    return 1ull << level;
  }

  auto hist_at_shard(uint64_t const shard) {
    DCHECK(hist_array::is_not_leveled);
    return hist_.get(shard, 0 /* TODO: 0 is arbitrary here*/);
  }
  auto hist_at_shard(uint64_t const shard) const {
    DCHECK(hist_array::is_not_leveled);
    return hist_.get(shard, 0 /* TODO: 0 is arbitrary here*/);
  }

  auto hist_at_level(uint64_t const level) {
    DCHECK(hist_array::is_not_sharded);
    return hist_.get(0 /* TODO: 0 is arbitrary here*/, level);
  }
  auto hist_at_level(uint64_t const level) const {
    DCHECK(hist_array::is_not_sharded);
    return hist_.get(0 /* TODO: 0 is arbitrary here*/, level);
  }

  auto borders_at_shard(uint64_t const shard) {
    DCHECK(borders_array::is_not_leveled);
    return borders_.get(shard, 0 /* TODO: 0 is arbitrary here*/);
  }
  auto borders_at_shard(uint64_t const shard) const {
    DCHECK(borders_array::is_not_leveled);
    return borders_.get(shard, 0 /* TODO: 0 is arbitrary here*/);
  }

  auto borders_at_level(uint64_t const level) {
    DCHECK(borders_array::is_not_sharded);
    return borders_.get(0 /* TODO: 0 is arbitrary here*/, level);
  }
  auto borders_at_level(uint64_t const level) const {
    DCHECK(borders_array::is_not_sharded);
    return borders_.get(0 /* TODO: 0 is arbitrary here*/, level);
  }

  auto hist_at(uint64_t const shard, uint64_t const level) {
    return hist_.get(shard, level);
  }
  auto hist_at(uint64_t const shard, uint64_t const level) const {
    return hist_.get(shard, level);
  }
  auto borders_at(uint64_t const shard, uint64_t const level) {
    return borders_.get(shard, level);
  }
  auto borders_at(uint64_t const shard, uint64_t const level) const {
    return borders_.get(shard, level);
  }

  uint64_t rho(size_t level, size_t i) {
    return rho_.get(level, i);
  }

  void set_rho(size_t level, size_t i, uint64_t val) {
    rho_.set(level, i, val);
  }

  static bool constexpr compute_zeros = !is_tree;
  static bool constexpr compute_rho = rho_type<is_tree>::compute_rho;

  // TODO: return span
  std::vector<uint64_t>& zeros() {
    return zeros_;
  }

  bv_type<require_initialization>& bv() {
    return bv_;
  }

  bv_type<require_initialization> const& bv() const {
    return bv_;
  }

  // TODO: Remove or make functional
  void discard_non_merge_data() {
    // Not used in merge algorithm
  }

private:
  hist_array hist_;
  borders_array borders_;
  std::vector<uint64_t> zeros_;
  bv_type<require_initialization> bv_;
  rho_type<is_tree> rho_;
  uint64_t levels_;
}; // ctx_generic

/******************************************************************************/
