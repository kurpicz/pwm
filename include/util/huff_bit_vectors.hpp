/*******************************************************************************
 * include/util/huff_bit_vectors.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "common.hpp"

class huff_bit_vectors {
public:
  huff_bit_vectors() data_(0), level_sizes_(0) { }

  huff_bit_vectors(const std::vector<uint64_t>& level_sizes) :
    data_(level_sizes.size()), level_sizes_(level_sizes) {

    uint64_t total_size = 0;
    for (const auto level_size : level_sizes_) {
      total_size += level_size;
    }

    data[0] = new uint64_t[word_size(total_size)];
    memset(data_[0], 0,
      (word_size(total_size) * sizeof(uint64_t)));

    for (uint64_t level = 1; level < level_sizes_.size(); ++level) {
      data_[level] = data_[level - 1] + word_size(level_sizes_[level - 1]);
    }
  }
  
  // non-copyable
  huff_bit_vectors(const huff_bit_vectors&) = delete;
  huff_bit_vectors& operator =(const huff_bit_vectors&) = delete;

  huff_bit_vectors(huff_bit_vectors&& other) : data_(std::move(other.data_)),
    level_sizes_(std::move(other.level_sizes_)) { }

  huff_bit_vectors& operator =(huff_bit_vectors&& other) {
    if (data_.size() > 0) {
      delete[] data_[0];
    }
    data_ = std::move(other.data_);
    level_sizes_ = std::move(other.level_sizes_);

    return *this;
  }

  ~huff_bit_vectors() {
    if  (data_.size() > 0) {
      delete[] data_[0];
    }
  }

  inline uint64_t levels() const {
    return data_.size();
  }

  inline std::vector<uint64_t> level_sizes() const {
    return level_sizes_;
  }

  inline uint64_t level_size(const uint64_t level) const {
    return level_sizes_[level];
  }

  inline const std::vector<uint64_t*>& vec() const {
    return data_;
  }

  inline std::vector<uint64_t*>& vec() {
    return data_;
  }

private:
  std::vector<uint64_t*> data_;
  std::vector<uint64_t> level_sizes_;
};

/******************************************************************************/
