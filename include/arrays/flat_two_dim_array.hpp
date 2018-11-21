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

template <typename IndexType>
class base_flat_two_dim_array {
public:
  base_flat_two_dim_array() = default;

  base_flat_two_dim_array(const uint64_t levels)
      : levels_(levels), data_(levels + 1), level_bit_sizes_(levels) {
    assert(levels > 0);
  }

  base_flat_two_dim_array(base_flat_two_dim_array&& other) = default;
  base_flat_two_dim_array& operator=(base_flat_two_dim_array&& other) {
    if (*this != other) {
      if (data_.size() > 0) {
        delete[] data_[0];
      }
      data_ = std::move(other.data_);
      levels_ = other.levels_;
      other.levels_ = 0;
      level_bit_sizes_ = std::move(other.level_bit_sizes_);
    }
    return *this;
  }

  ~base_flat_two_dim_array() {
    if (data_.size() > 0) {
      delete[] data_[0];
    }
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
            data_[0] == other.data_[0]);
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
    auto ptr = data_[index];
    auto nptr = data_[index + 1];
    return {ptr, size_t(nptr - ptr)};
  }

  inline span<IndexType> operator[](const uint64_t index) {
    DCHECK(index < levels());
    auto ptr = data_[index];
    auto nptr = data_[index + 1];
    return {ptr, size_t(nptr - ptr)};
  }

protected:
  inline const std::vector<IndexType*>& raw_data() const {
    return data_;
  }

  inline std::vector<IndexType*>& raw_data() {
    return data_;
  }

  uint64_t levels_ = 0;
  std::vector<IndexType*> data_;
  std::vector<uint64_t> level_bit_sizes_;

}; // class base_flat_two_dim_array

template <typename IndexType, class size_function,
          bool requires_initialization = true>
class flat_two_dim_array : public base_flat_two_dim_array<IndexType> {
  using base = base_flat_two_dim_array<IndexType>;

public:
  flat_two_dim_array() : base::base_flat_two_dim_array() {}

  template <typename... SizeFunctionArgs>
  flat_two_dim_array(const uint64_t levels, SizeFunctionArgs... size_f_args)
      : base::base_flat_two_dim_array(levels) {

    auto& level_bit_sizes_ = base::level_bit_sizes_;
    auto& data_ = base::data_;

    uint64_t data_size = 0;
    for (uint64_t level = 0; level < levels; ++level) {
      const uint64_t level_size =
          size_function::level_size(level, size_f_args...);
      // If its a bit vector, we still want to knwo how many bits there are
      // actually stored in each level, not just the number of computer words.
      if constexpr (size_function::is_bit_vector) {
        level_bit_sizes_[level] = level_size;
        data_size += word_size(level_size);
      } else {
        level_bit_sizes_[level] = level_size * (sizeof(IndexType) >> 3);
        data_size += level_size;
      }
    }
    data_[0] = new IndexType[data_size];

    if constexpr (requires_initialization) {
      #pragma omp parallel
      {
        const uint64_t omp_rank = omp_get_thread_num();
        const uint64_t omp_size = omp_get_num_threads();
        assert(omp_size == shards);

        const uint64_t local_size =
          (data_size / omp_size) + ((omp_rank < data_size % omp_size) ? 1 : 0);
        const uint64_t offset = (omp_rank * (data_size / omp_size)) +
                                 std::min<uint64_t>(omp_rank,
                                                    data_size % omp_size);
        memset(data_[0] + offset, 0, local_size * sizeof(IndexType));
      }
    }
    
    for (uint64_t level = 1; level < data_.size(); ++level) {
      const uint64_t level_size =
          size_function::level_size(level - 1, size_f_args...);
      if constexpr (size_function::is_bit_vector) {
        data_[level] = data_[level - 1] + word_size(level_size);
      } else {
        data_[level] = data_[level - 1] + level_size;
      }
    }
  }

}; // class flat_two_dim_array

/******************************************************************************/
