/*******************************************************************************
 * include/util/flat_two_dim_array.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstring>
#include <iostream>
#include "util/stxxl_helper.hpp"

template <typename IndexType, class size_function>
class flat_two_dim_array_external {

public:
  flat_two_dim_array_external() {}

  template <typename... SizeFunctionArgs>
  flat_two_dim_array_external(const uint64_t levels, SizeFunctionArgs... size_f_args)
  : levels_(levels), data_(levels + 1) {
    assert(levels > 0);
    level_sizes_.reserve(levels + 1);
    level_offsets_.reserve(levels + 1);
    uint64_t data_size = 0;
    for (uint64_t level = 0; level < levels; ++level) {
      uint64_t level_size = size_function::level_size(level, size_f_args...);
      level_sizes_.push_back(level_size);
      level_offsets_.push_back(data_size);
      data_size += level_size;
    }
    level_sizes_.push_back(0);
    level_offsets_.push_back(data_size);
    data_size_ = data_size;
    data_.resize(data_size);
  }

  flat_two_dim_array_external(flat_two_dim_array_external&& other)
  : levels_(other.levels_), data_size_(other.data_size_) { 
    level_sizes_.swap(other.level_sizes_);
    level_offsets_.swap(other.level_offsets_);
    zeros_.swap(other.zeros_);
    data_.swap(other.data_);
  }

  //~ flat_two_dim_array_external& operator =(flat_two_dim_array_external&& other) {
    //~ if (*this != other) {
      //~ if (data_.size() > 0) {
        //~ delete[] data_[0];
      //~ }
      //~ data_ = std::move(other.data_);
      //~ levels_ = other.levels_;
    //~ }
    //~ return *this;
  //~ }

  flat_two_dim_array_external(const flat_two_dim_array_external&) = delete;
  flat_two_dim_array_external& operator =(const flat_two_dim_array_external&) = delete;

  inline uint64_t levels() const {
    return levels_;
  }
  
  inline uint64_t data_size() const {
    return data_size_;
  }

  inline uint64_t level_size(uint64_t level) const {
    return level_sizes_[level];
  }

  inline std::vector<uint64_t> &level_sizes() {
    return level_sizes_;
  }

  inline std::vector<uint64_t> &level_offsets() {
    return level_offsets_;
  }

  inline const stxxlvector_offset<IndexType> operator [](const uint64_t level) const {
    return stxxlvector_offset<IndexType>(data_, level_offsets_[level]);
  }

  inline stxxlvector_offset<IndexType> operator [](const uint64_t level) {
    return stxxlvector_offset<IndexType>(data_, level_offsets_[level]);
  }

  inline const stxxlvector<IndexType>& raw_data() const {
    return data_;
  }

  inline stxxlvector<IndexType>& raw_data() {
    return data_;
  }
  
  inline std::vector<uint64_t> const& zeros() const {
    return zeros_;
  }
  
  inline void setZeros(const std::vector<uint64_t>& vec) {
    zeros_ = vec;
  }

private:
  uint64_t levels_ = 0;
  uint64_t data_size_ = 0;
  std::vector<uint64_t> level_sizes_;
  std::vector<uint64_t> level_offsets_;
  stxxlvector<IndexType> data_;
  
  std::vector<uint64_t> zeros_;
}; // class flat_two_dim_array

/******************************************************************************/
