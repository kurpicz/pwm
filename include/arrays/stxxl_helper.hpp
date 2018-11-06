/*******************************************************************************
 * include/util/common.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <stxxl/vector>
#include <stxxl/bits/io/linuxaio_file.h>
#include <stxxl/bits/io/syscall_file.h>
#include "util/type_for_bytes.hpp"

template <typename value_type>
using stxxlvector = typename stxxl::VECTOR_GENERATOR<value_type>::result;

template <typename value_type>
using stxxlreader = typename stxxl::VECTOR_GENERATOR<value_type>::result::bufreader_type;

template <typename value_type>
using stxxlwriter = typename stxxl::VECTOR_GENERATOR<value_type>::result::bufwriter_type;

template <int word_width>
using external_vector = stxxlvector<typename type_for_bytes<word_width>::type>;


class stxxl_files {
private:
  static constexpr int mode =
      stxxl::file::open_mode::RDWR |
      stxxl::file::open_mode::CREAT |
      stxxl::file::open_mode::TRUNC;

  std::vector<stxxl::syscall_file *> instance_files;
  std::vector<bool> instance_file_status;

  static stxxl_files& getInstance() {
    static stxxl_files instance;
    return instance;
  }

  static auto &files() {
    return getInstance().instance_files;
  }

  static auto &used() {
    return getInstance().instance_file_status;
  }

  stxxl_files() {};

  ~stxxl_files() {
    for(auto ptr : instance_files) {
      delete ptr;
    }
  }

public:
  stxxl_files(const stxxl_files &) = delete;
  void operator=(const stxxl_files &) = delete;

  static void addFile(const std::string &path) {
    files().push_back(new stxxl::syscall_file(path, mode));
    used().push_back(false);
  }

  // may return a non empty vector
  template <typename vector_type>
  static vector_type getVector(unsigned id, bool verbose = false) {
    bool exists = id < files().size();
    if(exists && !used()[id]) {
      return vector_type(files()[id]);
    } else {
      if(verbose) {
        std::cerr << "Trying to use EM buffer file with ID \"" << id << "\", ";
        if (!exists) std::cerr << "which has not been specified." << std::endl;
        else std::cerr << "which is already in use." << std::endl;
      }
      return vector_type();
    }
  }

  // only reset counter, if all vectors created by nextVector() have been destroyed
  static void reset_usage() {
    for(unsigned i = 0; i < used().size(); ++i)
      used()[i] = false;
  }
};

// simple accessor to simulate 2D array on 1D stxxlvector
// [for testing only, as [] operator is very expensive]
template <typename value_type>
class stxxlvector_offset {

private:
  const stxxlvector<value_type> &vec;
  const uint64_t off;

public:
  stxxlvector_offset(const stxxlvector<value_type> &vector, uint64_t offset)
  : vec(vector), off(offset) {}
  
  inline const value_type operator [](const uint64_t index) const {
    return vec[index + off];
  }
};

/******************************************************************************/
