/*******************************************************************************
 * include/util/external_flat_two_dim_array.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstring>
#include <iostream>

#include <stxxl/vector>

//~ #include <util/file_array.hpp>

template <typename IndexType, class size_function>
class external_flat_two_dim_array {
    
  using stxxlvector = typename stxxl::VECTOR_GENERATOR<IndexType>::result;
  //~ using filearray = file_array<IndexType>;
  //~ using defaultarray = std::vector<IndexType>;
  using em_array_type = stxxlvector;

public:
  external_flat_two_dim_array() { }

  template <typename... SizeFunctionArgs>
  external_flat_two_dim_array(const uint64_t levels, SizeFunctionArgs... size_f_args)
  : levels_(levels), level_offsets_(levels), level_sizes_(levels) {
    assert(levels > 0);
    
    uint64_t total_elements = 0;
    for (uint64_t level = 0; level < levels; ++level) {
      uint64_t current_elements = size_function::level_size(level, size_f_args...);
      level_offsets_[level] = total_elements;
      level_sizes_[level] = current_elements;
      total_elements += current_elements;
    }
    
    raw_data_.resize(total_elements);
  }

  external_flat_two_dim_array(external_flat_two_dim_array&& other)
  : levels_(other.levels_),
    raw_data_(std::forward<em_array_type>(other.raw_data_)),
    level_offsets_(std::forward<std::vector<uint64_t>>(other.level_offsets_)),
    level_sizes_(std::forward<std::vector<uint64_t>>(other.level_sizes_)) { }


  external_flat_two_dim_array& operator =(external_flat_two_dim_array&& other) {
    if (*this != other) {
      raw_data_ = std::move(other.raw_data_);
      level_offsets_ = std::move(other.level_offsets_);
      level_sizes_ = std::move(other.level_sizes_);
      levels_ = other.levels_;
    }
    return *this;
  }

  ~external_flat_two_dim_array() {
  }

  external_flat_two_dim_array(const external_flat_two_dim_array&) = delete;
  external_flat_two_dim_array& operator =(const external_flat_two_dim_array&) = delete;

  bool operator ==(const external_flat_two_dim_array& other) const {
    return 
      raw_data_ == other.raw_data_ &&
      level_offsets_ == other.level_offsets_ &&
      level_sizes_ == other.level_sizes_;
  }

  bool operator !=(const external_flat_two_dim_array& other) const {
    return !(*this == other);
  }

  inline uint64_t levels() const {
    return levels_;
  }

  inline uint64_t level_size(uint64_t level) const {
    return level_sizes_[level] * sizeof(IndexType);
  }

  inline decltype(em_array_type().cbegin()) operator [](const uint64_t index) const {
    return raw_data_.cbegin() + level_offsets_[index];
  }

  inline decltype(em_array_type().begin()) operator [](const uint64_t index) {
    return raw_data_.begin() + level_offsets_[index];
  }

  inline auto& raw_data() {
    return raw_data_;
  }
  
  inline auto& level_offsets() {
    return level_offsets_;
  }

private:

  uint64_t levels_;
  em_array_type raw_data_;
  std::vector<uint64_t> level_offsets_;
  std::vector<uint64_t> level_sizes_;
}; // class external_flat_two_dim_array

/******************************************************************************/
