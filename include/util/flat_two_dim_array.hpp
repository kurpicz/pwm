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

template <typename IndexType, class size_function>
class flat_two_dim_array {

public:
  flat_two_dim_array() : data_(0) { }

  template <typename... SizeFunctionArgs>
  flat_two_dim_array(const uint64_t levels, SizeFunctionArgs... size_f_args)
  : levels_(levels), data_(levels + 1), level_bit_sizes_(levels) {
    assert(levels > 0);
    uint64_t data_size = 0;
    for (uint64_t level = 0; level < levels; ++level) {
      const uint64_t level_size =
        size_function::level_size(level, size_f_args...);
      // If its a bit vector, we still want to knwo how many bits there are
      // actually stored in each level, not just
      if (size_function::is_bit_vector) { // TODO: C++17 if constexpr
        level_bit_sizes_[level] = level_size;
        data_size += word_size(level_size);
      } else {
        level_bit_sizes_[level] = level_size * (sizeof(IndexType) >> 3);
        data_size += level_size;
      }
    }
    data_[0] = new IndexType[data_size];
    memset(data_[0], 0, data_size * sizeof(IndexType));
    for (uint64_t level = 1; level < data_.size(); ++level) {
      const uint64_t level_size =
        size_function::level_size(level - 1, size_f_args...);
      if (size_function::is_bit_vector) { // TODO: C++17 if constexpr
        data_[level] = data_[level - 1] + word_size(level_size);
      } else {
        data_[level] = data_[level - 1] + level_size;
      }
    }
  }

  flat_two_dim_array(flat_two_dim_array&& other) = default;

  flat_two_dim_array& operator =(flat_two_dim_array&& other) {
    if (*this != other) {
      if (data_.size() > 0) {
        delete[] data_[0];
      }
      data_ = std::move(other.data_);
      levels_ = other.levels_;
    }
    return *this;
  }

  ~flat_two_dim_array() {
    if (data_.size() > 0) {
      delete[] data_[0];
    }
  }

  flat_two_dim_array(const flat_two_dim_array&) = delete;
  flat_two_dim_array& operator =(const flat_two_dim_array&) = delete;

  bool operator ==(const flat_two_dim_array& other) const {
    if ((data_.size() == 0 && other.data_.size() != 0) ||
        (data_.size() != 0 && other.data_.size() == 0)) {
      return false;
    }
    // Here we know that either both have length 0 or both have length > 0.
    return ((data_.size() == 0 && other.data_.size() == 0) ||
      data_[0] == other.data_[0]);
  }

  bool operator !=(const flat_two_dim_array& other) const {
    return !(*this == other);
  }

  inline uint64_t levels() const {
    return levels_;
  }

  inline uint64_t level_size(uint64_t level) const {
    return data_[level + 1] - data_[level];
  }

  inline uint64_t level_bit_size(uint64_t level) const {
    return level_bit_sizes_[level];
  }

  inline const IndexType* operator [](const uint64_t index) const {
    return data_[index];
  }

  inline IndexType* operator [](const uint64_t index) {
    return data_[index];
  }

  inline const std::vector<IndexType*>& raw_data() const {
    return data_;
  }

  inline std::vector<IndexType*>& raw_data() {
    return data_;
  }

private:
  uint64_t levels_;
  std::vector<IndexType*> data_;
  std::vector<uint64_t> level_bit_sizes_;

}; // class flat_two_dim_array

/******************************************************************************/
