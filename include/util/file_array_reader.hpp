/*******************************************************************************
 * include/util/file_stream.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <chrono>
#include <cstdlib>
#include <fcntl.h>
#include <tuple>
#include <unistd.h>
#include <errno.h>

#include "util/macros.hpp"

// meant for mostly sequential use only!
template <typename AlphabetType, uint64_t buffer_size = 10 * 1024 * 1024>
class file_array {

public:
  file_array(const std::string& file_name, const uint64_t size=0)
    : file_name_(file_name), buffer_first_pos_(0), size_(size) {

    buffer_ = new AlphabetType[buffer_size];

    // create new file
    file_descriptor_ = open(file_name.c_str(), O_RDWR | O_CREAT, S_IRWXU);
    if (file_descriptor_ == -1) std::cout << "Constructor FD" << std::endl; //std::exit(-1);
    
    // fill file with 0-bits
    if (size > 0) {
      if(ftruncate(file_descriptor_, (size + buffer_size) * sizeof(AlphabetType)) < 0)
        std::cout << "Constructor truncate" << std::endl; //std::exit(-1);
      //~ if(lseek(file_descriptor_, 0, SEEK_SET) < 0) std::cout << "Seek truncate" << std::endl; //std::exit(-1);
    }
  }
  
  file_array(const uint64_t size=0) 
    : file_array("output_tesla", 
    //std::string(
    //  std::chrono::duration_cast<std::chrono::milliseconds>(
    //    std::chrono::system_clock::now().time_since_epoch())),
      size) {}

  ~file_array() {
    write_buffer();
    //~ close(file_descriptor_);
  }
  
  void resize(uint64_t size) {
      size_ = size;
      if(lseek(file_descriptor_, 0, SEEK_SET) < 0) std::cout << "Resize seek1" << std::endl; //std::exit(-1);
      if(ftruncate(file_descriptor_, (size + buffer_size) * sizeof(AlphabetType)) < 0)
        std::cout << "Resize truncate" << std::endl; //std::exit(-1);
      //~ if(lseek(file_descriptor_, 0, SEEK_SET) < 0) std::cout << "Resize seek2 truncate" << std::endl; //std::exit(-1);
  }

  inline void write_buffer() {
    signed code;
    if ((code = lseek(file_descriptor_, buffer_first_pos_ * sizeof(AlphabetType), SEEK_SET)) < 0)
      std::cout << "Write seek" << errno << std::endl; //std::exit(-1);
    int64_t bytes_written = write(file_descriptor_, buffer_,
      buffer_size * sizeof(AlphabetType));
    if (PWM_UNLIKELY(bytes_written) == -1) std::cout << "Write bytes" << bytes_written << " " << buffer_first_pos_ << " " << size_ << " " << errno << std::endl; //std::exit(-1);
  }
  
  inline void read_buffer(const uint64_t index) {
    write_buffer();
    if (lseek(file_descriptor_, index * sizeof(AlphabetType), SEEK_SET) < 0) {
      std::cout << "Read seek" << std::endl; //std::exit(-1);
    }
    int64_t bytes_read = read(file_descriptor_, buffer_, buffer_size * sizeof(AlphabetType));
    if(bytes_read != buffer_size * sizeof(AlphabetType)) {
      std::cout << "Read bytes " << bytes_read << " " << buffer_size * sizeof(AlphabetType) << " " << index << std::endl;
      std::exit(-1);
    }
    buffer_first_pos_ = index;
  }

  std::string file_name() const {
    return file_name_;
  }

  AlphabetType& operator [](const uint64_t index) {
    if (PWM_UNLIKELY(index < buffer_first_pos_ ||
      index >= buffer_first_pos_ + buffer_size)) {
      read_buffer(index);
    }
    return buffer_[index % buffer_size];
  }
  
  const AlphabetType& operator [](const uint64_t index) const {
    if (PWM_UNLIKELY(index < buffer_first_pos_ ||
      index >= buffer_first_pos_ + buffer_size)) {
      read_buffer(index);
    }
    return buffer_[index % buffer_size];
  }

  bool operator ==(const file_array& other) const {
    return std::tie(buffer_, file_descriptor_) == std::tie(other.buffer_,
      other.file_descriptor_);
  }

  bool operator !=(const file_array& other) const {
    return *this != other;
  }
  
  file_array_iterator<AlphabetType, buffer_size> begin() {
    return file_array_iterator<AlphabetType, buffer_size>(*this, 0);
  }
  
  file_array_const_iterator<AlphabetType, buffer_size> begin() const {
    return file_array_const_iterator<AlphabetType, buffer_size>(*this, 0);
  }
  
  file_array_const_iterator<AlphabetType, buffer_size> cbegin() const {
    return file_array_const_iterator<AlphabetType, buffer_size>(*this, 0);
  }

private:
  std::string file_name_;
  uint64_t buffer_first_pos_;
  uint64_t size_;
  AlphabetType * buffer_;
  int32_t file_descriptor_;

}; // class file_array

/******************************************************************************/
