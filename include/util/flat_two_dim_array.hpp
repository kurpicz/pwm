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

template <typename IndexType, class size_function>
class flat_two_dim_array {

public:
  flat_two_dim_array() : data_(0) { }

  template <typename... SizeFunctionArgs>
  flat_two_dim_array(const uint64_t levels, SizeFunctionArgs... size_f_args)
  : data_(levels) {
    assert(levels > 0);
    uint64_t data_size = 0;
    for (uint64_t level = 0; level < levels; ++level) {
      data_size += size_function::level_size(level, size_f_args...);
    }
    data_[0] = new IndexType[data_size];
    memset(data_[0], 0, data_size * sizeof(IndexType));
    for (uint64_t level = 1; level < levels; ++level) {
      data_[level] = data_[level - 1] +
        size_function::level_size(level - 1, size_f_args...);
    }
  }

  flat_two_dim_array(flat_two_dim_array&& other)
  : data_(std::forward<std::vector<IndexType*>>(other.data_)) { }

  flat_two_dim_array& operator =(flat_two_dim_array&& other) {
    if (*this != other) {
      if (data_.size() > 0) {
        delete[] data_[0];
      }
      data_ = std::move(other.data_);
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
    return data_.size();
  }

  inline const IndexType* operator [](const uint64_t index) const {
    return data_[index];
  }

  inline IndexType* operator [](const uint64_t index) {
    return data_[index];
  }

  // TODO: Rename to data_vector
  inline const std::vector<IndexType*>& vec() const {
    return data_;
  }

  // TODO: Rename to data_vector
  inline std::vector<IndexType*>& vec() {
    return data_;
  }

private:
  std::vector<IndexType*> data_;

}; // class flat_two_dim_array

/******************************************************************************/