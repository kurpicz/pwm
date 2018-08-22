/*******************************************************************************
 * include/util/flat_two_dim_array.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstring>
#include <iostream>

#include "util/common.hpp"

template <typename IndexType>
class base_flat_two_dim_array {
public:
  base_flat_two_dim_array() : data_(0) { }

  base_flat_two_dim_array(const uint64_t levels)
  : levels_(levels), data_(levels + 1), level_bit_sizes_(levels) {
    assert(levels > 0);
  }

  base_flat_two_dim_array(base_flat_two_dim_array&& other) = default;

  base_flat_two_dim_array& operator =(base_flat_two_dim_array&& other) {
    if (*this != other) {
      if (data_.size() > 0) {
        delete[] data_.at(0);
      }
      data_ = std::move(other.data_);
      levels_ = other.levels_;
      level_bit_sizes_ = std::move(other.level_bit_sizes_);
    }
    return *this;
  }

  ~base_flat_two_dim_array() {
    if (data_.size() > 0) {
      delete[] data_.at(0);
    }
  }

  base_flat_two_dim_array(const base_flat_two_dim_array&) = delete;
  base_flat_two_dim_array& operator =(const base_flat_two_dim_array&) = delete;

  bool operator ==(const base_flat_two_dim_array& other) const {
    if ((data_.size() == 0 && other.data_.size() != 0) ||
        (data_.size() != 0 && other.data_.size() == 0)) {
      return false;
    }
    // Here we know that either both have length 0 or both have length > 0.
    return ((data_.size() == 0 && other.data_.size() == 0) ||
      data_.at(0) == other.data_.at(0));
  }

  bool operator !=(const base_flat_two_dim_array& other) const {
    return !(*this == other);
  }

  inline uint64_t levels() const {
    return levels_;
  }

  inline uint64_t level_size(uint64_t level) const {
    return data_.at(level + 1) - data_.at(level);
  }

  inline uint64_t level_bit_size(uint64_t level) const {
    return level_bit_sizes_.at(level);
  }

  inline const IndexType* operator [](const uint64_t index) const {
    return data_.at(index);
  }

  inline IndexType* operator [](const uint64_t index) {
    return data_.at(index);
  }

  inline const std::vector<IndexType*>& raw_data() const {
    return data_;
  }

  inline std::vector<IndexType*>& raw_data() {
    return data_;
  }

protected:
  uint64_t levels_;
  std::vector<IndexType*> data_;
  std::vector<uint64_t> level_bit_sizes_;

}; // class base_flat_two_dim_array

template <typename IndexType, class size_function>
class flat_two_dim_array: public base_flat_two_dim_array<IndexType> {
  using base = base_flat_two_dim_array<IndexType>;
public:
  flat_two_dim_array(): base::base_flat_two_dim_array() {}

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
      // actually stored in each level, not just
      if (size_function::is_bit_vector) { // TODO: C++17 if constexpr
        level_bit_sizes_.at(level) = level_size;
        data_size += word_size(level_size);
      } else {
        level_bit_sizes_.at(level) = level_size * (sizeof(IndexType) >> 3);
        data_size += level_size;
      }
    }
    data_.at(0) = new IndexType[data_size];
    memset(data_.at(0), 0, data_size * sizeof(IndexType));
    for (uint64_t level = 1; level < data_.size(); ++level) {
      const uint64_t level_size =
        size_function::level_size(level - 1, size_f_args...);
      if (size_function::is_bit_vector) { // TODO: C++17 if constexpr
        data_.at(level) = data_.at(level - 1) + word_size(level_size);
      } else {
        data_.at(level) = data_.at(level - 1) + level_size;
      }
    }
  }

}; // class flat_two_dim_array

/******************************************************************************/
