/*******************************************************************************
 * include/algorithm.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
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

// NB: These headers provide definitions used by the .cpp implementation files
// and the CONSTRUCTION_REGISTER macro. We flag them as being provided by this
// header to omit them from warning output
#include <string> // IWYU pragma: export
#include <cstdint> // IWYU pragma: export

#include "arrays/memory_types.hpp"
#include "arrays/stxxl_helper.hpp"
#include "construction/wavelet_structure.hpp"
#include "util/common.hpp"
#include "util/type_for_bytes.hpp"

template <typename in_type, typename out_type>
class construction_algorithm;

template <typename Algorithm>
class concrete_algorithm;

template <typename in_type, typename out_type>
class algorithm_list {
public:
  algorithm_list(const algorithm_list& other) = delete;
  algorithm_list(algorithm_list&& other) = delete;
  void operator=(const algorithm_list& other) = delete;
  void operator=(algorithm_list&& other) = delete;

  inline static algorithm_list& get_algorithm_list() {
    static algorithm_list list;
    return list;
  }

  inline void
  register_algorithm(construction_algorithm<in_type, out_type> const* algo) {
    algorithms_.push_back(algo);
  }

  inline auto begin() {
    return algorithms_.cbegin();
  }
  inline auto end() {
    return algorithms_.cend();
  }

  template <typename filter_function>
  inline auto filtered(filter_function f) {
    std::vector<construction_algorithm<in_type, out_type> const*> r;
    for (auto e : *this) {
      if (f(e)) {
        r.push_back(e);
      }
    }
    return r;
  }

private:
  algorithm_list() {}

  // List of static pointers to the different algorithms
  std::vector<construction_algorithm<in_type, out_type> const*> algorithms_;
}; // class algorithm_list

template <typename in_type, typename out_type>
class construction_algorithm {
public:
  construction_algorithm(std::string const& name, std::string const& description)
      : name_(name), description_(description) {
    algorithm_list<in_type, out_type>::get_algorithm_list().register_algorithm(
        this);
  }

  virtual out_type compute_bitvector(const in_type& global_text,
                                     const uint64_t size,
                                     const uint64_t levels) const = 0;

  inline float median_time(const in_type& global_text,
                           const uint64_t size,
                           const uint64_t levels,
                           const uint64_t runs) const {
    std::vector<float> times;
    for (uint64_t run = 0; run < runs; ++run) {
      auto begin_time = std::chrono::high_resolution_clock::now();
      compute_bitvector(global_text, size, levels);
      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = end_time - begin_time;
      auto millis =
          std::chrono::duration_cast<std::chrono::milliseconds>(duration);
      times.emplace_back(static_cast<float>(millis.count()));
    }
    std::sort(times.begin(), times.end());
    return (times[(runs - 1) >> 1] + times[(runs >> 1)]) / 2;
  }

  inline void memory_peak(const in_type& global_text,
                          const uint64_t size,
                          const uint64_t levels) const {
    compute_bitvector(global_text, size, levels);
  }

  virtual bool is_parallel() const = 0;
  virtual bool is_tree() const = 0;
  virtual uint8_t word_width() const = 0;
  virtual bool is_huffman_shaped() const = 0;
  virtual uint64_t huffman_bit_size(const in_type& global_text,
                                    const uint64_t size,
                                    const uint64_t levels) const = 0;
  virtual bool is_input_external() const = 0;
  virtual bool is_output_external() const = 0;

  std::string name() const {
    return name_;
  }

  std::string description() const {
    return description_;
  }

  inline void print_info() const {
    std::cout << name_ << ": " << description_ << std::endl;
  }

private:
  std::string name_;
  std::string description_;
}; // class construction_algorithm

// TODO: merge concrete algorithm, construction algorithm

template <typename Algorithm>
class concrete_algorithm
    : construction_algorithm<typename algo_type<Algorithm>::in,
                             typename algo_type<Algorithm>::out> {
public:
  using in_type = typename algo_type<Algorithm>::in;
  using out_type = typename algo_type<Algorithm>::out;

  concrete_algorithm(std::string const& name, std::string const& description)
      : construction_algorithm<in_type, out_type>(name, description) {}

  out_type compute_bitvector(const in_type& global_text,
                             const uint64_t size,
                             const uint64_t levels) const override {
    return Algorithm::compute(global_text, size, levels);
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

  inline uint64_t huffman_bit_size([[maybe_unused]] const in_type& global_text,
                                   [[maybe_unused]] const uint64_t size,
                                   [[maybe_unused]] const uint64_t levels)
    const override {

    uint64_t bit_size = 0;
    if constexpr (Algorithm::is_huffman_shaped) {
      auto&& wx = compute_bitvector(global_text, size, levels);
      auto&& bit_vecs = wx.bvs();
      for (uint64_t i = 0; i < bit_vecs.levels(); ++i) {
        bit_size += bit_vecs.level_bit_size(i);
      }
    }
    return bit_size;
  }

  bool is_input_external() const override {
    return Algorithm::external_in;
  }

  bool is_output_external() const override {
    return Algorithm::external_out;
  }

}; // class concrete_algorithm

#define CONSTRUCTION_REGISTER(algo_name, algo_description, ws)                 \
  static const auto _cstr_algo_##ws##_register =                               \
      concrete_algorithm<ws>(algo_name, algo_description);

#endif // ALGORITHM_HEADER

/******************************************************************************/
