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

  inline static algorithm_list& get_algorithm_list() {
    static algorithm_list list;
    return list;
  }

  inline void register_algorithm(construction_algorithm const* algo) {
    algorithms_.push_back(algo);
  }

  inline auto begin() { return algorithms_.cbegin(); }
  inline auto end() { return algorithms_.cend(); }

private:
  algorithm_list() { }

  // List of static pointers to the different algorithms
  std::vector<construction_algorithm const*> algorithms_;
}; // class algorithm_list

class construction_algorithm {
public:
  construction_algorithm(std::string name, std::string description)
    : name_(name), description_(description) {
    algorithm_list::get_algorithm_list().register_algorithm(this);
  }

  inline virtual std::pair<std::vector<uint64_t*>, std::vector<uint64_t>> compute_bitvector(
        const void* global_text, const uint64_t size, const uint64_t levels) const = 0;
  inline virtual bool is_parallel() const  = 0;
  inline virtual bool is_tree() const  = 0;
  inline virtual uint8_t word_width() const  = 0;


  inline void print_info() const {
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

  inline std::pair<std::vector<uint64_t*>, std::vector<uint64_t>> compute_bitvector(
    const void* global_text, const uint64_t size, const uint64_t levels) const override
  {
    using text_vec_type =
      std::vector<typename type_for_bytes<WaveletStructure::word_width>::type>;

    auto const* text = static_cast<text_vec_type const*>(global_text);
    ws_ = WaveletStructure(*text, size, levels);
    return ws_.get_bv_and_zeros();
  }

  bool is_parallel() const override {
    return WaveletStructure::is_parallel;
  }

  bool is_tree() const override {
    return WaveletStructure::is_tree;
  }

  uint8_t word_width() const override {
    return WaveletStructure::word_width;
  }

private:
  WaveletStructure ws_;

}; // class concrete_algorithm

#define CONSTRUCTION_REGISTER(algo_name, algo_description, wavelet_structure) \
  static const auto _cstr_algo_ ## wavelet_structure ## _register            \
    = concrete_algorithm<wavelet_structure>(algo_name, algo_description);

#endif // ALGORITHM_HEADER

/******************************************************************************/
