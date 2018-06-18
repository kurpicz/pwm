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

#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <vector>

#include "util/common.hpp"
#include "util/file_stream.hpp"
#include "util/type_for_bytes.hpp"
#include "util/wavelet_structure.hpp"

class construction_algorithm_base;

template<bool is_semi_external>
class construction_algorithm;

template<typename Algorithm, bool is_semi_external>
class concrete_algorithm;

template<bool is_semi_external>
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

  inline void register_algorithm(construction_algorithm<is_semi_external> const* algo) {
    algorithms_.push_back(algo);
  }

  inline auto begin() { return algorithms_.cbegin(); }
  inline auto end() { return algorithms_.cend(); }

private:
  algorithm_list() { }

  // List of static pointers to the different algorithms
  std::vector<construction_algorithm<is_semi_external> const*> algorithms_;
}; // class algorithm_list


class construction_algorithm_base {
public:

  construction_algorithm_base(std::string name, std::string description)
    : name_(name), description_(description) {}

  virtual float median_time(const void* global_text, const uint64_t size,
      const uint64_t levels, size_t runs) const = 0;
  virtual void memory_peak(const void* global_text, const uint64_t size,
      const uint64_t levels ) const = 0;
  virtual bool is_parallel() const = 0;
  virtual bool is_tree() const = 0;
  virtual uint8_t word_width() const = 0;
  virtual bool is_huffman_shaped() const = 0;

  std::string name() const {
    return name_;
  }

  std::string description() const {
    return description_;
  }

  inline void print_info() const {
    std::cout << name_ << ": " << description_ << std::endl;
  }

protected:
  std::string name_;
  std::string description_;
}; // class construction_algorithm

template<>
class construction_algorithm<false> : public construction_algorithm_base {
public:
  construction_algorithm(std::string name, std::string description)
    : construction_algorithm_base(name, description) {
    algorithm_list<false>::get_algorithm_list().register_algorithm(this);
  }

  virtual wavelet_structure<false> compute_bitvector(const void* global_text,
    const uint64_t size, const uint64_t levels) const = 0;
}; // class construction_algorithm

template<>
class construction_algorithm<true> : public construction_algorithm_base {
public:
  construction_algorithm(std::string name, std::string description)
    : construction_algorithm_base(name, description) {
    algorithm_list<true>::get_algorithm_list().register_algorithm(this);
  }

  virtual wavelet_structure<true> compute_bitvector(const void* global_text,
    const uint64_t size, const uint64_t levels) const = 0;
}; // class construction_algorithm


template <typename Algorithm, bool is_semi_external>
class concrete_algorithm : public construction_algorithm<is_semi_external> {

public:
  concrete_algorithm(std::string name, std::string description)
  : construction_algorithm<is_semi_external>(name, description) { }

  inline wavelet_structure<is_semi_external> compute_bitvector(const void* global_text,
    const uint64_t size, const uint64_t levels) const override {
    using text_vec_type =
      std::vector<typename type_for_bytes<Algorithm::word_width>::type>;
    auto const* text = static_cast<text_vec_type const*>(global_text);
    return Algorithm::compute(text->data(), size, levels);
  }

  float median_time(const void* global_text, const uint64_t size,
      const uint64_t levels, size_t runs) const override {
    std::vector<float> times;
    using text_vec_type =
      std::vector<typename type_for_bytes<Algorithm::word_width>::type>;
    const auto* text = static_cast<const text_vec_type*>(global_text);
    for (size_t run = 0; run < runs; ++run) {
      auto begin_time = std::chrono::high_resolution_clock::now();
      Algorithm::compute(text->data(), size, levels);
      auto end_time = std::chrono::high_resolution_clock::now();
      times.emplace_back(static_cast<float>(
        std::chrono::duration_cast<std::chrono::milliseconds>(
          end_time - begin_time).count()));
    }
    std::sort(times.begin(), times.end());
    return times[runs >> 1];
  }

  void memory_peak(const void* global_text, const uint64_t size,
    const uint64_t levels) const override {
    using text_vec_type =
      std::vector<typename type_for_bytes<Algorithm::word_width>::type>;
    const auto* text = static_cast<const text_vec_type*>(global_text);
    Algorithm::compute(text->data(), size, levels);
  }

  bool is_parallel() const override {
    return Algorithm::is_parallel;
  }

  bool is_tree() const override {
    return Algorithm::is_tree;
  }

  uint8_t word_width() const override {
    return Algorithm::word_width;
  }

  bool is_huffman_shaped() const override {
    return Algorithm::is_huffman_shaped;
  }

}; // class concrete_algorithm

#define CONSTRUCTION_REGISTER(algo_name, algo_description, ws) \
  static const auto _cstr_algo_ ## ws ## _register             \
    = concrete_algorithm<ws, false>(algo_name, algo_description);
    
#define CONSTRUCTION_REGISTER_SE(algo_name, algo_description, ws) \
  static const auto _cstr_algo_ ## ws ## _register             \
    = concrete_algorithm<ws, true>(algo_name, algo_description);

#endif // ALGORITHM_HEADER

/******************************************************************************/
