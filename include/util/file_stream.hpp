/*******************************************************************************
 * include/util/file_stream.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstdlib>
#include <fcntl.h>
#include <unistd.h>

#include "util/macros.hpp"

template <typename AlphabetType, uint64_t buffer_size=1024*1024>
class file_stream {

public:
  file_stream(const std::string& file_name, const uint64_t offset=0)
  : offset_(offset), max_text_pos_(0) {

    file_descriptor_ = open(file_name.c_str(), O_RDONLY);
    if (file_descriptor_ == -1) {
      std::exit(-1);
    }
    // We will only scan the text in 
    posix_fadvise(file_descriptor_, 0, 0, 1);
    reset_stream();
    read_next_buffer();
    --max_text_pos_; // Position starts at 0; we add buffer_size elements.
  }

  ~file_stream() {
    close(file_descriptor_);
  }

  void reset_stream() {
    if (lseek(file_descriptor_, offset_, SEEK_SET) < 0) {
      std::exit(-1);
    }
  }

  const AlphabetType operator [](const uint64_t index) {
    if (PWM_UNLIKELY(index > max_text_pos_)) {
      read_next_buffer();
    }
    return buffer_[index % buffer_size];
  }

private:
  inline void read_next_buffer() {
    int64_t bytes_read = read(file_descriptor_, buffer_.data(),
      buffer_size * sizeof(AlphabetType));
    if (PWM_UNLIKELY(bytes_read) == -1) {
      std::exit(-1);
    }
    max_text_pos_ += bytes_read / sizeof(AlphabetType);
  }

private:
  uint64_t offset_;
  uint64_t max_text_pos_;
  std::array<AlphabetType, buffer_size> buffer_;
  int32_t file_descriptor_;
}; // class file_stream

/******************************************************************************/
