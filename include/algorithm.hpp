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

#include "wavelet_structures.hpp"

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

  virtual void run() = 0;
  virtual bool is_tree() = 0;

  void print_info() {
    std::cout << name_ << ": " << description_ << std::endl;
  }

private:
  std::string name_;
  std::string description_;
}; // class construction_algorithm

template <typename AlphabetType,
          typename SizeType,
          template <typename> class WaveletStructure>
class concrete_construction_algorithm : construction_algorithm {

public:
  using construct_function = void (*)(const AlphabetType* text,
    const SizeType size, const SizeType levels);

  concrete_construction_algorithm(std::string name, std::string description,
    construct_function cstr_fnct) : construction_algorithm(name, description),
    cstr_fnct_(cstr_fnct) { }

  inline void run() { }

  bool is_tree() {
    return WaveletStructure<SizeType>::is_tree;
  }

private:
  construct_function cstr_fnct_;
};

#define CONSTRUCTION_ALGORITHM_REGISTER(algo_name, algo_description, \
  a_t, s_t, w_t, cstr_fct)                                           \
  static const auto* _cstr_algo_ ## cstr_fct ## _register            \
    = new concrete_construction_algorithm<a_t, s_t, w_t>(            \
      algo_name, algo_description, cstr_fct);

#endif // ALGORITHM_HEADER

/******************************************************************************/
