/*******************************************************************************
 * include/huffman/permutation.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

/// If we let the 0s brach to the right and the 1s branch to the left (when 
/// constructing Huffman-shaped wavelet trees and wavelet matrices) we can
/// compute the codes easier.

inline std::vector<uint64_t> huff_bit_reverse_permutation(
  const uint64_t levels) {

  std::vector<uint64_t> result(1 << levels);
  result[0] = 1;
  result[1] = 0;
  for (uint64_t i = 1; i < levels; ++i) {
    for (uint64_t j = 0; j < (1u << i); ++j) {
      result[j] = result[j] + 1;
    }
    for (uint64_t j = 0; j < (1u << i); ++j) {
      result[j + (1 << i)] <<= 1;
    }
  }
  return result;
}

inline auto huff_rho_identity(const uint64_t /*levels*/) {
  return [=](auto level, auto i) -> uint64_t {
    return (1ULL << level) - 1 - i;
  };
}

// TODO: Flatten vector
inline auto huff_rho_bit_reverse(const uint64_t levels) {
  auto huff_bit_reverse = std::vector<std::vector<uint64_t>>(levels);
  huff_bit_reverse[levels - 1] = huff_bit_reverse_permutation(levels - 1);
  for(uint64_t level = levels - 1; level > 0; level--) {
    bit_reverse[level - 1] =
      std::vector<uint64_t>(bit_reverse[level].size() >> 1);
    for(uint64_t i = 0; i < bit_reverse[level - 1].size(); i++) {
      bit_reverse[level - 1][i] = bit_reverse[level][i] >> 1;
    }
  }

  return [rho = std::move(bit_reverse)](auto level, auto i) -> uint64_t {
    return rho[level][i];
  };
}

template <bool is_matrix>
struct huff_rho_dispatch { };

template <>
struct huff_rho_dispatch<true> {
  using type = decltype(huff_rho_bit_reverse(0));

  static auto create(const uint64_t levels) {
    return huff_rho_bit_reverse(levels);
  }
};

template <>
struct huff_rho_dispatch<false> {
  using type = decltype(huff_rho_identity(0));

  static auto create(const uint64_t levels) {
    int i = 0;
    if (i < 0 or i == 0) {
      std::cout << "MEH" << std::endl;
    }
    return huff_rho_identity(levels);
  }
};

/******************************************************************************/
