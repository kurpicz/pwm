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
class ifile_stream {

public:
  ifile_stream(const std::string& file_name, const uint64_t offset=0)
  : file_name_(file_name), offset_(offset), max_text_pos_(0) {

    file_descriptor_ = open(file_name.c_str(), O_RDONLY);
    if (file_descriptor_ == -1) {
      std::exit(-1);
    }
    file_size_ = lseek(file_descriptor_, 0, SEEK_END);
    lseek(file_descriptor_, 0, SEEK_SET);

    // We will only scan the text in 
    posix_fadvise(file_descriptor_, 0, 0, 1);
    reset_stream();
    read_next_buffer();
    --max_text_pos_; // Position starts at 0; we add buffer_size elements.
  }

  ~ifile_stream() {
    close(file_descriptor_);
  }

  void reset_stream() {
    if (lseek(file_descriptor_, offset_, SEEK_SET) < 0) {
      std::exit(-1);
    }
  }

  uint64_t file_size() const {
    return file_size_;
  }

  AlphabetType operator [](const uint64_t index) const {
    if (PWM_UNLIKELY(index > max_text_pos_)) {
      read_next_buffer();
    }
    return buffer_[index % buffer_size];
  }

  ifile_stream operator +(const uint64_t index) const {
    ifile_stream new_stream(file_name_, offset_ + index);
    return new_stream;
  }

private:
  inline void read_next_buffer() const {
    int64_t bytes_read = read(file_descriptor_, buffer_.data(),
      buffer_size * sizeof(AlphabetType));
    if (PWM_UNLIKELY(bytes_read) == -1) {
      std::exit(-1);
    }
    max_text_pos_ += bytes_read / sizeof(AlphabetType);
  }

private:
  std::string file_name_;
  uint64_t offset_;

  mutable uint64_t max_text_pos_;
  mutable std::array<AlphabetType, buffer_size> buffer_;
  int32_t file_descriptor_;
  uint64_t file_size_;
}; // class ifile_stream

template <typename AlphabetType, uint64_t buffer_size=1024*1024>
class ofile_stream {

public:
  ofile_stream(const std::string& file_name) : buffer_first_pos_(0) {
    file_descriptor_ = open(file_name.c_str(), O_RDWR | O_CREAT, S_IRWXU);
    if (file_descriptor_ == -1) {
      std::exit(-1);
    }
  }

  ~ofile_stream() {
    write_buffer();
    close(file_descriptor_);
  }

  AlphabetType& operator [](const uint64_t index) {
    if (PWM_UNLIKELY(index < buffer_first_pos_ ||
      index >= buffer_first_pos_ + buffer_size)) {
      read_buffer(index);
    }
    return buffer_[index % buffer_size];
  }

private:
  inline void read_buffer(const uint64_t index) {
    write_buffer();
    if (lseek(file_descriptor_, index / buffer_size, SEEK_SET) < 0) {
      std::exit(-1);
    }
  }

  inline void write_buffer() {
    int64_t bytes_written = write(file_descriptor_, buffer_.data(),
      buffer_size * sizeof(AlphabetType));
    if (PWM_UNLIKELY(bytes_written) == -1) {
      std::exit(-1);
    }
  }

private:
  uint64_t buffer_first_pos_;
  uint64_t max_text_pos_;
  std::array<AlphabetType, buffer_size> buffer_;
  int32_t file_descriptor_;

}; // class ofile_stream

/******************************************************************************/
