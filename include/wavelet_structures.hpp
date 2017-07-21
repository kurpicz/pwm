#pragma once
#ifndef WAVELET_STRUCTURES_HEADER
#define WAVELET_STRUCTURES_HEADER


#include "common.hpp"

template <permutation_type<> Permutation, bool IsTree>
struct wavelet_structure {
  static constexpr bool is_tree = IsTree;
  static constexpr permutation_type<> permutation = Permutation;
}; // struct wavelet_structure

using wavelet_tree = wavelet_structure<identity_function, true>;
using wavelet_matrix = wavelet_structure<bit_reverse_permutation, false>;

#endif // WAVELET_STRUCTURES_HEADER
