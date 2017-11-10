/*******************************************************************************
 * include/util/bit_vectors.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "common.hpp"

class bit_vectors {

public:
  bit_vectors() : size_(0) { }

  bit_vectors(const uint64_t size, const uint64_t levels) : data_(levels),
    size_(size) {
    assert(levels != 0);

    data_[0] = new uint64_t[word_size(size) * levels];
    memset(data_[0], 0, (word_size(size) * sizeof(uint64_t)) * levels);

    for (uint64_t level = 1; level < levels; ++level) {
      data_[level] = data_[level - 1] + word_size(size);
    }
  }

  //non-copyable
  bit_vectors(const bit_vectors&) = delete;
  bit_vectors& operator =(const bit_vectors&) = delete;

  bit_vectors(bit_vectors&& other) : data_(std::move(other.data_)),
    size_(other.size_) { }

  bit_vectors& operator =(bit_vectors&& other) {
    if (data_.size() > 0) {
      delete[] data_[0];
    }
    data_ = std::move(other.data_);
    size_ = other.size_;

    return *this;
  }

  ~bit_vectors() {
    if (data_.size() > 0) {
      delete[] data_[0];
    }
  }

  inline uint64_t levels() const {
    return data_.size();
  }

  inline uint64_t size() const {
    return size_;
  }

  inline const std::vector<uint64_t*>& vec() const {
    return data_;
  }

  inline std::vector<uint64_t*>& vec() {
    return data_;
  }

private:
  std::vector<uint64_t*> data_;
  uint64_t size_;
};

/******************************************************************************/
