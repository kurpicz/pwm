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

#include "util/common.hpp"
#include "util/file_stream.hpp"
#include "util/type_for_bytes.hpp"
#include "util/memory_types.hpp"
#include "util/wavelet_structure.hpp"
#include "util/stxxl_helper.hpp"

template<memory_mode mem_mode>
class construction_algorithm;

template<typename Algorithm>
class concrete_algorithm;

template<memory_mode mem_mode>
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

  inline void register_algorithm(construction_algorithm<mem_mode> const* algo) {
    algorithms_.push_back(algo);
  }

  inline auto begin() { return algorithms_.cbegin(); }
  inline auto end() { return algorithms_.cend(); }

private:
  algorithm_list() { }

  // List of static pointers to the different algorithms
  std::vector<construction_algorithm<mem_mode> const*> algorithms_;
}; // class algorithm_list


template<memory_mode mem_mode>
class construction_algorithm {
protected:
  std::string name_;
  std::string description_;

  static constexpr bool is_input_external_ = 
      mem_mode == memory_mode::external || 
      mem_mode == memory_mode::external_input;
  
  static constexpr bool is_output_external_ =
      mem_mode == memory_mode::external || 
      mem_mode == memory_mode::external_output;
      
public:

  construction_algorithm(std::string name, std::string description)
    : name_(name), description_(description) {
    algorithm_list<mem_mode>::get_algorithm_list().register_algorithm(this);
  }

  virtual float median_time(const void* global_text, const uint64_t size,
      const uint64_t levels, size_t runs) const = 0;
  virtual void memory_peak(const void* global_text, const uint64_t size,
      const uint64_t levels ) const = 0;
  virtual bool is_parallel() const = 0;
  virtual bool is_tree() const = 0;
  virtual uint8_t word_width() const = 0;
  virtual bool is_huffman_shaped() const = 0;
  
  constexpr bool is_input_external() const {
    return 
      mem_mode == memory_mode::external || 
      mem_mode == memory_mode::external_input;
  }
  
  constexpr bool is_output_external() const {
    return 
      mem_mode == memory_mode::external || 
      mem_mode == memory_mode::external_output;
  }
  
  memory_mode mmode() {
    return mem_mode;
  }

  std::string name() const {
    return name_;
  }

  std::string description() const {
    return description_;
  }

  inline void print_info() const {
    std::cout << name_ << ": " << description_ << std::endl;
  }
  
  virtual typename out_type<is_output_external_>::type compute_bitvector(const void * global_text,
    const uint64_t size, const uint64_t levels) const = 0;

}; // class construction_algorithm



template <typename Algorithm>
class concrete_algorithm : public construction_algorithm<Algorithm::mem_mode>{
private:
  
  using input_type = 
    typename in_type<construction_algorithm<
      Algorithm::mem_mode>::is_input_external_, 
      Algorithm::word_width
    >::type;
    
  using output_type = 
    typename out_type<construction_algorithm<
      Algorithm::mem_mode>::is_output_external_
    >::type;
  
public:
  concrete_algorithm(std::string name, std::string description)
  : construction_algorithm<Algorithm::mem_mode>(name, description) { }

  inline output_type compute_bitvector(const void* global_text,
    const uint64_t size, const uint64_t levels) const override {
    auto const &input = * static_cast<input_type const *>(global_text);
    return Algorithm::template compute<input_type>(input, size, levels);
  }

  float median_time(const void* global_text, const uint64_t size,
      const uint64_t levels, size_t runs) const override {
    std::vector<float> times;
    const input_type &input = *static_cast<input_type const *>(global_text);
    for (size_t run = 0; run < runs; ++run) {
      auto begin_time = std::chrono::high_resolution_clock::now();
      Algorithm::template compute<input_type>(input, size, levels);
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
    compute_bitvector(global_text, size, levels);
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
    = concrete_algorithm<ws>(algo_name, algo_description);

#endif // ALGORITHM_HEADER

/******************************************************************************/
