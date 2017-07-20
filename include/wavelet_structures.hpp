#pragma once
#ifndef WAVELET_STRUCTURES_HEADER
#define WAVELET_STRUCTURES_HEADER

#include <type_traits>

#include "common.hpp"

template <typename SizeType,
          permutation_type<SizeType> Permutation,
          bool tree>
class wavelet_structure {
public:

  template <bool _is_tree = tree>
  wavelet_structure(const std::enable_if_t<_is_tree, SizeType> levels)
  : bit_vector(levels), zeros(0) { }

  template <bool _is_tree = tree>
  wavelet_structure(const std::enable_if_t<!_is_tree, SizeType> levels)
  : bit_vector(levels), zeros(levels, 0) { }

  static constexpr bool is_tree = tree;
  static constexpr permutation_type<SizeType> permutation = Permutation;

public:
  std::vector<uint64_t*> bit_vector;
  std::vector<SizeType>  zeros;
};

template <typename SizeType>
using wavelet_tree = wavelet_structure<SizeType, identity_function, true>;

template <typename SizeType>
using wavelet_matrix = wavelet_structure<SizeType, bit_reverse_permutation, false>;

#endif // WAVELET_STRUCTURES_HEADER
