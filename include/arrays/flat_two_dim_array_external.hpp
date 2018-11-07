/*******************************************************************************
 * include/util/flat_two_dim_array.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "stxxl_helper.hpp"
#include <cstring>
#include <iostream>

template <typename IndexType, class size_function>
class flat_two_dim_array_external {

public:
  flat_two_dim_array_external() {}

  flat_two_dim_array_external(const uint64_t levels,
                              const uint64_t size,
                              const int fileid = -1)
      : levels_(levels) {
    assert(levels > 0);
    if (fileid > -1)
      data_ = stxxl_files::getVectorPermanent<stxxlvector<IndexType>>(fileid);
    data_.clear();
    level_sizes_.reserve(levels + 1);
    level_offsets_.reserve(levels + 1);
    uint64_t data_size = 0;
    const unsigned bits_per_element = 8 * sizeof(IndexType);
    for (uint64_t level = 0; level < levels; ++level) {
      uint64_t level_size = size_function::level_size(level, size);
      level_sizes_.push_back(level_size);
      level_offsets_.push_back(data_size);
      data_size += (level_size + bits_per_element - 1) / bits_per_element;
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

  flat_two_dim_array_external(const flat_two_dim_array_external&) = delete;
  flat_two_dim_array_external&
  operator=(const flat_two_dim_array_external&) = delete;

  inline uint64_t levels() const {
    return levels_;
  }

  inline uint64_t data_size() const {
    return data_size_;
  }

  inline uint64_t level_size(uint64_t level) const {
    return level_sizes_[level];
  }

  inline std::vector<uint64_t>& level_sizes() {
    return level_sizes_;
  }

  inline std::vector<uint64_t>& level_offsets() {
    return level_offsets_;
  }

  inline const stxxlvector_offset<IndexType>
  operator[](const uint64_t level) const {
    return stxxlvector_offset<IndexType>(data_, level_offsets_[level]);
  }

  inline stxxlvector_offset<IndexType> operator[](const uint64_t level) {
    return stxxlvector_offset<IndexType>(data_, level_offsets_[level]);
  }

  /*  inline const stxxlvector<IndexType>& raw_data() const {
      return data_;
    }*/

  inline stxxlvector<IndexType>& raw_data() {
    return data_;
  }

  inline void set_raw_data(stxxlvector<IndexType>& raw) {
    data_ = raw;
  }

  inline std::vector<uint64_t>& zeros() {
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
