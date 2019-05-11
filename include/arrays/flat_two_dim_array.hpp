/*******************************************************************************
 * include/arrays/flat_two_dim_array.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstring>
#include <iostream>
#include <omp.h>

#include "span.hpp"
#include "util/common.hpp"
#include "util/debug_assert.hpp"

template <typename IndexType>
class base_flat_two_dim_array {
public:
  base_flat_two_dim_array() = default;

  base_flat_two_dim_array(const uint64_t levels)
      : levels_(levels), data_(levels), level_bit_sizes_(levels) {
    DCHECK(levels > 0);
  }

  base_flat_two_dim_array(base_flat_two_dim_array&& other) {
    move_from(std::move(other));
  }
  base_flat_two_dim_array& operator=(base_flat_two_dim_array&& other) {
    move_from(std::move(other));
    return *this;
  }

  ~base_flat_two_dim_array() {
    destroy();
  }

  base_flat_two_dim_array(const base_flat_two_dim_array&) = delete;
  base_flat_two_dim_array& operator=(const base_flat_two_dim_array&) = delete;

  bool operator==(const base_flat_two_dim_array& other) const {
    if ((data_.size() == 0 && other.data_.size() != 0) ||
        (data_.size() != 0 && other.data_.size() == 0)) {
      return false;
    }
    // Here we know that either both have length 0 or both have length > 0.
    return ((data_.size() == 0 && other.data_.size() == 0) ||
            data_[0].data() == other.data_[0].data());
  }

  bool operator!=(const base_flat_two_dim_array& other) const {
    return !(*this == other);
  }

  inline uint64_t levels() const {
    return levels_;
  }

  inline uint64_t level_bit_size(uint64_t level) const {
    return level_bit_sizes_[level];
  }

  inline span<IndexType const> operator[](const uint64_t index) const {
    DCHECK(index < levels());
    return data_[index];
  }

  inline span<IndexType> operator[](const uint64_t index) {
    DCHECK(index < levels());
    return data_[index];
  }

protected:
  uint64_t levels_ = 0;
  std::vector<span<IndexType>> data_;
  std::vector<uint64_t> level_bit_sizes_;
  IndexType* allocation_ = nullptr;

  inline void destroy() {
    if (data_.size() > 0) {
      DCHECK(allocation_ != nullptr);
      delete[] allocation_;
    }
  }

  inline void move_from(base_flat_two_dim_array&& other) {
    destroy();

    data_ = std::move(other.data_);
    level_bit_sizes_ = std::move(other.level_bit_sizes_);

    levels_ = other.levels_;
    other.levels_ = 0;

    allocation_ = other.allocation_;
    other.allocation_ = nullptr;
  }

}; // class base_flat_two_dim_array

/// Aa 2D vector that stores the the different levels in a continous
/// allocation.
///
/// Each level might be seprated by some padding to ensure they
/// start in different cache lines.
template <typename IndexType, class config>
class flat_two_dim_array : public base_flat_two_dim_array<IndexType> {
  using base = base_flat_two_dim_array<IndexType>;

  static_assert((CACHELINE_SIZE % sizeof(IndexType)) == 0);

  inline static uint64_t
  compute_elements_till_next_cacheline(uint64_t initial_idx) {
    auto remainder = initial_idx % elements_per_cacheline();
    if (remainder == 0) {
      remainder = elements_per_cacheline();
    }
    return elements_per_cacheline() - remainder;
  }
  inline static uint64_t elements_per_cacheline() {
    return CACHELINE_SIZE / sizeof(IndexType);
  }
  struct level_data_size_info {
    uint64_t level_array_size;
    uint64_t level_bit_size;
    uint64_t level_array_size_plus_cacheline_align;
  };
  template <typename... SizeFunctionArgs>
  inline static level_data_size_info
  compute_level_size(uint64_t level, SizeFunctionArgs... size_f_args) {
    const uint64_t level_size = config::level_size(level, size_f_args...);
    level_data_size_info ret;

    // If its a bit vector, we still want to knwo how many bits there are
    // actually stored in each level, not just the number of computer words.
    if constexpr (config::is_bit_vector) {
      ret.level_bit_size = level_size;
      ret.level_array_size = word_size(level_size);
    } else {
      ret.level_bit_size = level_size * (sizeof(IndexType) >> 3);
      ret.level_array_size = level_size;
    }

    auto extra = compute_elements_till_next_cacheline(ret.level_array_size);

    ret.level_array_size_plus_cacheline_align = ret.level_array_size + extra;

    return ret;
  }

public:
  flat_two_dim_array() : base::base_flat_two_dim_array() {}

  template <typename... SizeFunctionArgs>
  flat_two_dim_array(const uint64_t levels, SizeFunctionArgs... size_f_args)
      : base::base_flat_two_dim_array(levels) {

    auto& level_bit_sizes_ = base::level_bit_sizes_;
    auto& data_ = base::data_;
    DCHECK(data_.size() == levels);

    // compute element counts and allocation size per level
    uint64_t data_size = 0;
    for (uint64_t level = 0; level < levels; ++level) {
      const auto info = compute_level_size(level, size_f_args...);

      level_bit_sizes_[level] = info.level_bit_size;
      data_size += info.level_array_size_plus_cacheline_align;
    }

    // We need to align the allocation to cachelines,
    // so add some extra wiggle room
    data_size += elements_per_cacheline();

    // allocate and optionally initialize the data
    auto& allocation_ = base::allocation_;
    allocation_ = new IndexType[data_size];
    if constexpr (config::requires_initialization) {
      #pragma omp parallel
      {
        const uint64_t omp_rank = omp_get_thread_num();
        const uint64_t omp_size = omp_get_num_threads();

        const uint64_t local_size = (data_size / omp_size) +
                                    ((omp_rank < data_size % omp_size) ? 1 : 0);
        const uint64_t offset =
            (omp_rank * (data_size / omp_size)) +
            std::min<uint64_t>(omp_rank, data_size % omp_size);
        memset(allocation_ + offset, 0, local_size * sizeof(IndexType));
      }
    }

    // align first level pointer to the start of the next cacheline
    auto aligned_ptr = allocation_;
    auto ptr_addr = uint64_t(allocation_);
    DCHECK(ptr_addr % sizeof(IndexType) == 0);
    aligned_ptr +=
        compute_elements_till_next_cacheline(ptr_addr / sizeof(IndexType));

    // setup the spans pointing to the individual levels in the allocation
    for (uint64_t level = 0; level < levels; ++level) {
      const auto info = compute_level_size(level, size_f_args...);

      data_.at(level) = span<IndexType>(aligned_ptr, info.level_array_size);
      aligned_ptr += info.level_array_size_plus_cacheline_align;
    }
  }

}; // class flat_two_dim_array

/******************************************************************************/
