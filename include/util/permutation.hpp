/*******************************************************************************
 * include/util/permutation.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstdint>

using permutation_type = std::vector<uint64_t> (*)(const uint64_t levels);

inline std::vector<uint64_t> bit_reverse_permutation(const uint64_t levels) {
  std::vector<uint64_t> result(1 << levels);
  result[0] = 0;
  result[1] = 1;
  for (uint64_t i = 1; i < levels; ++i) {
    for (uint64_t j = 0; j < (1u << i); ++j) {
      result[j] <<= 1;
    }
    for (uint64_t j = 0; j < (1u << i); ++j) {
      result[j + (1 << i)] = result[j] + 1;
    }
  }
  return result;
}

inline auto rho_identity(const uint64_t /*levels*/) {
  return [](const auto /*level*/, const auto i) -> uint64_t { return i; };
}

// TODO: Flatten vector
inline auto rho_bit_reverse(const uint64_t levels) {
  auto bit_reverse = std::vector<std::vector<uint64_t>>(levels);
  bit_reverse[levels - 1] = bit_reverse_permutation(levels - 1);
  for (uint64_t level = levels - 1; level > 0; level--) {
    bit_reverse[level - 1] =
        std::vector<uint64_t>(bit_reverse[level].size() >> 1);
    for (uint64_t i = 0; i < bit_reverse[level - 1].size(); i++) {
      bit_reverse[level - 1][i] = bit_reverse[level][i] >> 1;
    }
  }

  return [rho = std::move(bit_reverse)](const auto level,
                                        const auto i) -> uint64_t {
    return rho[level][i];
  };
}

template <bool is_tree>
struct rho_dispatch {};

template <>
struct rho_dispatch<true> {
  using type = decltype(rho_identity(0));

  static auto create(uint64_t levels) {
    return rho_identity(levels);
  }
};

template <>
struct rho_dispatch<false> {
  using type = decltype(rho_bit_reverse(0));

  static auto create(uint64_t levels) {
    return rho_bit_reverse(levels);
  }
};

/******************************************************************************/
