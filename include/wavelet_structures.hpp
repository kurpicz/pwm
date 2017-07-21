#pragma once
#ifndef WAVELET_STRUCTURES_HEADER
#define WAVELET_STRUCTURES_HEADER

#include <type_traits>

#include "common.hpp"

template <permutation_type<> Permutation, bool IsTree>
class wavelet_structure {
public:
  static constexpr bool is_tree = IsTree;
  static constexpr permutation_type<> permutation = Permutation;
};

using wavelet_tree = wavelet_structure<identity_function, true>;
using wavelet_matrix = wavelet_structure<bit_reverse_permutation, false>;

#endif // WAVELET_STRUCTURES_HEADER
