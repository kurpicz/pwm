/*******************************************************************************
 * include/algorithm.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef ALGORITHM_HEADER
#define ALGORITHM_HEADER

#include <iostream>
#include <memory>
#include <vector>

#include "util/common.hpp"
#include "util/type_for_bytes.hpp"

class construction_algorithm;

class algorithm_list {
public:
  algorithm_list(const algorithm_list& other) = delete;
  algorithm_list(algorithm_list&& other) = delete;
  void operator = (const algorithm_list& other) = delete;
  void operator = (algorithm_list&& other) = delete;

  static algorithm_list& get_algorithm_list() {
    static algorithm_list list;
    return list;
  }

  void register_algorithm(std::unique_ptr<construction_algorithm>&& algo) {
    algorithms_.emplace_back(std::move(algo));
  }

  auto begin() { return algorithms_.begin(); }
  auto end() { return algorithms_.end(); }

private:
  algorithm_list() { }
  std::vector<std::unique_ptr<construction_algorithm>> algorithms_;
}; // class algorithm_list

class construction_algorithm {
public:
  construction_algorithm(std::string name, std::string description)
    : name_(name), description_(description) {
    algorithm_list::get_algorithm_list().register_algorithm(
      std::move(std::unique_ptr<construction_algorithm>(this)));
  }

  virtual std::pair<Bvs, std::vector<uint64_t>> compute_bitvector(
    const void* global_text, const uint64_t size, const uint64_t levels) = 0;
  virtual bool is_parallel() = 0;
  virtual bool is_tree() = 0;
  virtual uint8_t word_width() = 0;


  void print_info() {
    std::cout << name_ << ": " << description_ << std::endl;
  }

private:
  std::string name_;
  std::string description_;
}; // class construction_algorithm

template <typename WaveletStructure>
class concrete_algorithm : construction_algorithm {

public:
  concrete_algorithm(std::string name, std::string description)
  : construction_algorithm(name, description) { }

  inline auto compute_bitvector(const void* global_text, const uint64_t size,
    const uint64_t levels) {
    using text_type =
      typename type_for_bytes<WaveletStructure::word_width>::type;
    std::vector<text_type> text = *(static_cast<text_type*>(global_text));
    ws_ = WaveletStructure(text, size, levels);
    return ws_.get_bv_and_zeros();
  }

  bool is_parallel() {
    return WaveletStructure::is_parallel;
  }

  bool is_tree() {
    return WaveletStructure::is_tree;
  }

  uint8_t word_width() {
    return WaveletStructure::word_width;
  }

private:
  WaveletStructure ws_;

}; // class concrete_algorithm

#define CONSTRUCTION_REGISTER(algo_name, algo_description, wavelet_structure) \
  static const auto* _cstr_algo_ ## wavelet_structure ## _register            \
    = new concrete_algorithm<wavelet_structure>(algo_name, algo_description);

#endif // ALGORITHM_HEADER

/******************************************************************************/
