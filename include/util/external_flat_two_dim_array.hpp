/*******************************************************************************
 * include/util/external_flat_two_dim_array.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstring>
#include <iostream>

#include "util/file_stream.hpp"

template <typename IndexType, class size_function>
class external_flat_two_dim_array {

public:
  external_flat_two_dim_array() : data_(0) { }

  template <typename... SizeFunctionArgs>
  external_flat_two_dim_array(const uint64_t levels,
    SizeFunctionArgs... size_f_args) : levels_(levels), start_pos_(levels + 1) {

    assert(levels > 0);
    uint64_t data_size = 0;
    for (uint64_t level = 0; level < levels; ++level) {
      data_size += size_function::level_size(level, size_f_args...);
    }
    // data_[0] = new IndexType[data_size];
    // memset(data_[0], 0, data_size * sizeof(IndexType));
    start_pos_[0] = 0;
    for (uint64_t level = 1; level < start_pos_.size(); ++level) {
      start_pos_[level] = start_pos_[level - 1] +
        size_function::level_size(level - 1, size_f_args...);
    }
  }

  external_flat_two_dim_array(external_flat_two_dim_array&& other)
  : levels_(other.levels_),
    data_(std::forward<std::vector<IndexType*>>(other.data_)) { }

  external_flat_two_dim_array& operator =(external_flat_two_dim_array&& other) {
    if (*this != other) {
      if (data_.size() > 0) {
        delete[] data_[0];
      }
      data_ = std::move(other.data_);
      levels_ = other.levels_;
    }
    return *this;
  }

  ~external_flat_two_dim_array() {
    if (data_.size() > 0) {
      delete[] data_[0];
    }
  }

  external_flat_two_dim_array(const external_flat_two_dim_array&) = delete;
  external_flat_two_dim_array& operator =(
    const external_flat_two_dim_array&) = delete;

  bool operator ==(const external_flat_two_dim_array& other) const {
    return data_ == other.data_;
  }

  bool operator !=(const external_flat_two_dim_array& other) const {
    return !(*this == other);
  }

  inline uint64_t levels() const {
    return levels_;
  }

  inline uint64_t level_size(uint64_t level) const {
    return start_pos_[level + 1] - start_pos_[level];
  }

  inline const IndexType* operator [](const uint64_t index) const {
    return data_[index];
  }

  inline IndexType* operator [](const uint64_t index) {
    return data_[index];
  }

  inline const ofile_stream<IndexType>& raw_data() const {
    return data_;
  }

  inline ofile_stream<IndexType>& raw_data() {
    return data_;
  }

private:
  uint64_t levels_;
  std::vector<uint64_t> start_pos_;
  ofile_stream<IndexType> data_;

}; // class external_flat_two_dim_array

/******************************************************************************/
