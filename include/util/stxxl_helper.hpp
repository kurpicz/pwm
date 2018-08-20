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


//~ class stxxlvector_factory {
//~ public:
  //~ template <typename value_type>
  //~ static stxxlvector<value_type> createOnDisk() {
    //~ // TODO: Implement factory method to create vector on specific device
    //~ return stxxlvector<value_type>()
  //~ }
//~ };


class stxxl_files {
private:
  static constexpr int mode =
      stxxl::file::open_mode::RDWR |
      stxxl::file::open_mode::CREAT |
      stxxl::file::open_mode::TRUNC;

  unsigned instance_next_id = 0;
  std::vector<stxxl::syscall_file *> instance_files;

  static stxxl_files& getInstance() {
    static stxxl_files instance;
    return instance;
  }

  static auto &next_id() {
    return getInstance().instance_next_id;
  }

  static auto &files() {
    return getInstance().instance_files;
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
  }

  template <typename vector_type>
  static vector_type nextVector() {
    if(next_id() < files().size()) {
      return vector_type(files()[next_id()++]);
    } else {
      return vector_type();
    }
  }

  // only reset counter, if all vectors created by nextVector() have been destroyed
  static void resetCounter() {
    next_id() = 0;
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
  stxxlvector_offset(stxxlvector<value_type> &vector, uint64_t offset) 
  : vec(vector), off(offset) {}
  
  inline const value_type operator [](const uint64_t index) const {
    return vec[index + off];
  }
};

/******************************************************************************/
