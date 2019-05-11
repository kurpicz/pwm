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
#include "util/stats.hpp"

template <bool ext_in, bool ext_out, int word_width>
class construction_algorithm;

template <typename Algorithm>
class concrete_algorithm;

template <bool ext_in, bool ext_out, int word_width>
class algorithm_list {
  using calgo_type = construction_algorithm<ext_in, ext_out, word_width>;
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
  register_algorithm(calgo_type const* algo) {
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
    std::vector<calgo_type const*> r;
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
  std::vector<calgo_type const*> algorithms_;
}; // class algorithm_list

template <bool ext_in, bool ext_out, int width>
class construction_algorithm {
  using in_type = typename input_type<ext_in, width>::type;
  using out_type = typename output_type<ext_out>::type;
  using out_type_internal = typename output_type<false>::type;

  using stats_type = statistics<ext_in | ext_out>;
public:
  construction_algorithm(std::string const& name, std::string const& description)
      : name_(name), description_(description) {
    algorithm_list<ext_in, ext_out, width>::get_algorithm_list().
        register_algorithm(this);
  }

  virtual out_type compute_bitvector(const in_type& global_text,
                                     const uint64_t size,
                                     const uint64_t levels,
                                     stats_type& stats) const = 0;

  out_type compute_bitvector(const in_type& global_text,
                             const uint64_t size,
                             const uint64_t levels) const {
    stats_type dummy;
    return compute_bitvector(global_text, size, levels, dummy);
  }

  out_type_internal compute_internal_bitvector(const in_type& global_text,
                                               const uint64_t size,
                                               const uint64_t levels) const {
    if constexpr (ext_out)
      return compute_bitvector(global_text, size, levels)
          .getInternalStructure();
    else
      return compute_bitvector(global_text, size, levels);
  }

  inline stats_type median_time_stats(const in_type& global_text,
                                      const uint64_t size,
                                      const uint64_t levels,
                                      const uint64_t runs) const {
    std::vector<stats_type> stats(runs);
    for (uint64_t run = 0; run < runs; ++run) {
      stats[run].start();
      compute_bitvector(global_text, size, levels, stats[run]);
      stats[run].finish();
    }
    std::sort(stats.begin(), stats.end());
    return stats[(runs - 1) >> 1];
  }

  inline void memory_peak(const in_type& global_text,
                          const uint64_t size,
                          const uint64_t levels) const {
    compute_bitvector(global_text, size, levels);
  }

  virtual bool is_parallel() const = 0;
  virtual bool is_tree() const = 0;
  virtual bool is_huffman_shaped() const = 0;

  constexpr bool is_input_external() const { return ext_in; }
  constexpr bool is_output_external() const { return ext_out; }
  constexpr uint8_t word_width() const {return width; }

  virtual uint64_t huffman_bit_size(const in_type& global_text,
                                    const uint64_t size,
                                    const uint64_t levels) const = 0;


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
    : construction_algorithm<
        Algorithm::external_in,
        Algorithm::external_out,
        Algorithm::word_width> {
public:
  using in_type = typename input_type<Algorithm::external_in,
                                      Algorithm::word_width>::type;
  using out_type = typename output_type<Algorithm::external_out>::type;

  using base_class = construction_algorithm<
      Algorithm::external_in,
      Algorithm::external_out,
      Algorithm::word_width>;

  using stats_type = statistics<Algorithm::external_in | Algorithm::external_out>;

  concrete_algorithm(std::string name, std::string description)
      : base_class(name, description) {}

  out_type compute_bitvector(const in_type& global_text,
                             const uint64_t size,
                             const uint64_t levels,
                             stats_type &stats) const {
    if constexpr (Algorithm::stats_support)
      return Algorithm::compute(global_text, size, levels, stats);
    else
      return Algorithm::compute(global_text, size, levels);
  }

  bool is_parallel() const override {
    return Algorithm::is_parallel;
  }

  bool is_tree() const override {
    return Algorithm::is_tree;
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
      auto&& wx = base_class::compute_bitvector(global_text, size, levels);
      auto&& bit_vecs = wx.bvs();
      for (uint64_t i = 0; i < bit_vecs.levels(); ++i) {
        bit_size += bit_vecs.level_bit_size(i);
      }
    }
    return bit_size;
  }
}; // class concrete_algorithm

#define CONSTRUCTION_REGISTER(algo_name, algo_description, ws)                 \
  static const auto _cstr_algo_##ws##_register =                               \
      concrete_algorithm<ws>(algo_name, algo_description);

#define CONSTRUCTION_REGISTER_MEMBER(algo_name, algo_description, ws)          \
  const concrete_algorithm<ws> _cstr_algo_##ws##_register =                    \
      concrete_algorithm<ws>(algo_name, algo_description);

#endif // ALGORITHM_HEADER

/******************************************************************************/
