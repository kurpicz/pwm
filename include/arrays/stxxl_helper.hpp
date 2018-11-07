/*******************************************************************************
 * include/util/common.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <chrono>
#include <sstream>
#include <stxxl/bits/io/linuxaio_file.h>
#include <stxxl/bits/io/syscall_file.h>
#include <stxxl/vector>
#include "util/type_for_bytes.hpp"
#include "util/filesystem_util.hpp"

template <typename value_type>
using stxxlvector = typename stxxl::VECTOR_GENERATOR<value_type>::result;

template <typename value_type>
using stxxlreader =
    typename stxxl::VECTOR_GENERATOR<value_type>::result::bufreader_type;

template <typename value_type>
using stxxlwriter =
    typename stxxl::VECTOR_GENERATOR<value_type>::result::bufwriter_type;

template <int word_width>
using external_vector = stxxlvector<typename type_for_bytes<word_width>::type>;

class stxxl_files {
private:
  static constexpr int mode = stxxl::file::open_mode::RDWR |
                              stxxl::file::open_mode::CREAT |
                              stxxl::file::open_mode::TRUNC;

  std::string instance_id;
  std::vector<std::string> instance_directories;
  std::vector<stxxl::syscall_file*> instance_files;

  static stxxl_files& getInstance() {
    static stxxl_files instance;
    return instance;
  }

  static auto& files() {
    return getInstance().instance_files;
  }

  static auto& directories() {
    return getInstance().instance_directories;
  }

  static auto& id() {
    return getInstance().instance_id;
  }

  static std::string toBase36(uint64_t id) {
    std::ostringstream result;
    while(id > 0) {
      char c = id % 36;
      result << ((c < 10) ? char(c + 48) : char(c - 10 + 97));
      id /= 36;
    }
    return result.str();
  }

  stxxl_files(){
    auto now = std::chrono::high_resolution_clock::now();
    auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);
    auto epoch = now_ms.time_since_epoch();
    instance_id = toBase36(epoch.count());
  };

  ~stxxl_files() {
    for (auto ptr : instance_files) {
      delete ptr;
    }
  }

public:
  stxxl_files(const stxxl_files&) = delete;
  void operator=(const stxxl_files&) = delete;

  static bool addDirectory(const std::string& path, bool verbose = false) {
    directories().push_back(path);
    if(!createDirectory(path)) {
      if(verbose) std::cerr << "Cannot use given EM directory: " << path << std::endl;
      return false;
    }
    else if(verbose) std::cout << "Using EM directory: " << path << std::endl;
    return true;
  }

  // may return a non empty vector
  template <typename vector_type>
  static vector_type getVector(unsigned dir_id, bool verbose = false) {
    if(dir_id < directories().size()) {
      std::ostringstream file_stream;
      file_stream << directories()[dir_id] << "/" << id() << "."
                  << "vec" << sizeof(typename vector_type::value_type) * 8 << "."
                  << files().size();
      std::string file = file_stream.str();
      stxxl::syscall_file* stxxl_file = new stxxl::syscall_file(file, mode);
      files().push_back(stxxl_file);
      return vector_type(stxxl_file);
    } else {
      if (verbose) {
        std::cerr << "Trying to use EM directory with ID \"" << dir_id << "\", ";
        std::cerr << "which has not been specified." << std::endl;
      }
      return vector_type();
    }
  }
};

// simple accessor to simulate 2D array on 1D stxxlvector
// [for testing only, as [] operator is very expensive]
template <typename value_type>
class stxxlvector_offset {

private:
  const stxxlvector<value_type>& vec;
  const uint64_t off;

public:
  stxxlvector_offset(const stxxlvector<value_type>& vector, uint64_t offset)
      : vec(vector), off(offset) {}

  inline const value_type operator[](const uint64_t index) const {
    return vec[index + off];
  }
};

/******************************************************************************/
