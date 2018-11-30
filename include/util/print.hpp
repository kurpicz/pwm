/*******************************************************************************
 * include/util/debug.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <functional>
#include <iostream>
#include <string>

#include "construction/wavelet_structure.hpp"
#include "util/debug.hpp"
#include "util/debug_assert.hpp"

[[maybe_unused]] static void print_structure(std::ostream& out,
                                             const base_bit_vectors& bv,
                                             const std::vector<uint64_t>& zeros,
                                             bool is_tree,
                                             bool gap_words = false) {
  for (uint64_t i = 0; i < bv.levels(); i++) {
    out << "   bv[" << i << "]";
    out << "[";
    for (uint64_t j = 0; j < bv.level_bit_size(i); j++) {
      if (gap_words && j != 0 && (j % 64) == 0) {
        out << " ";
      }
      out << uint64_t(bit_at(bv[i], j)) << "";
    }
    out << "]";
    if (!is_tree && i != (bv.levels() - 1)) {
      out << " zeros[" << i << "] = " << zeros.at(i);
    }
    out << std::endl;
  }
  out << std::endl;
}

[[maybe_unused]] static void print_structure(std::ostream& out,
                                             wavelet_structure const& structure,
                                             bool gap_words = false) {
  const base_bit_vectors& bv = structure.bvs();
  const std::vector<uint64_t>& zeros = structure.zeros();
  print_structure(out, bv, zeros, structure.is_tree(), gap_words);
}

[[maybe_unused]] static auto get_block_split(size_t i,
                                             const base_bit_vectors& bv,
                                             uint64_t left,
                                             uint64_t right) {
  struct BlockSplit {
    uint64_t left;
    uint64_t mid;
    uint64_t right;
  };

  if (left > bv.level_bit_size(i)) {
    left = bv.level_bit_size(i);
  }
  DCHECK(left <= bv.level_bit_size(i));

  if (right > bv.level_bit_size(i)) {
    right = bv.level_bit_size(i);
  }
  DCHECK(right <= bv.level_bit_size(i));

  DCHECK(left <= right);

  uint64_t ones = 0;
  for (uint64_t j = left; j < right; j++) {
    bool bit = bit_at(bv[i], j);
    ones += bit;
  }
  uint64_t zeros = (right - left) - ones;

  return BlockSplit{left, left + zeros, right};
}

[[maybe_unused]] static void
print_tree(std::ostream& out, const base_bit_vectors& bv, bool padded = false) {
  struct Block {
    std::vector<bool> bits;
    uint64_t left;
    uint64_t right;
  };
  std::vector<std::vector<Block>> bvs(bv.levels());

  auto copy_bits = [&bv](auto i, auto& layer) {
    auto l = layer.left;
    auto r = layer.right;

    DCHECK(l <= r);
    auto v = std::vector<bool>();

    for (uint64_t j = l; j < r; j++) {
      if (j < bv.level_bit_size(i)) {
        v.push_back(bit_at(bv[i], j));
      }
    }

    layer.bits = std::move(v);
  };

  bvs.at(0).push_back({});
  bvs.at(0).at(0).left = 0;
  bvs.at(0).at(0).right = bv.level_bit_size(0);
  copy_bits(0, bvs.at(0).at(0));

  for (uint64_t i = 1; i < bv.levels(); i++) {
    auto& layer = bvs.at(i);

    for (auto& pblock : bvs.at(i - 1)) {
      auto res = get_block_split(i - 1, bv, pblock.left, pblock.right);

      auto lr = [&](auto l, auto r) {
        DCHECK(l <= r);
        layer.push_back({
            {},
            l,
            r,
        });
        copy_bits(i, layer.back());
      };

      lr(res.left, res.mid);
      lr(res.mid, res.right);
    }
  }
  for (uint64_t i = 0; i < bv.levels(); i++) {
    auto& layer = bvs.at(i);
    for (auto& block : layer) {
      auto pad = [&] {
        if (padded)
          for (size_t k = 0; k < (1ull << (bv.levels() - i - 1)); k++) {
            out << " ";
          }
      };
      pad();
      out << "[";
      for (bool bit : block.bits) {
        out << int(bit);
      }
      out << "]";
      pad();
    }
    out << std::endl;
  }
  out << std::endl;
}

template <typename T>
struct force_integer_trait {
  inline static void print(std::ostream& out, T const& v) {
    out << v;
  }
};
template <>
struct force_integer_trait<uint8_t> {
  inline static void print(std::ostream& out, uint8_t const& v) {
    out << uint64_t(v & 0xff);
  }
};
/// Print T, but convert it to int first if it is a char-like type
struct print_force_type {
  template <typename T>
  inline void operator()(std::ostream& out, T const& t) const {
    force_integer_trait<T>::print(out, t);
  }
};

template <typename list_type, typename map_type = print_force_type>
static std::ostream& print_list(std::ostream& out,
                                list_type const& list,
                                bool nice = false,
                                map_type fmt = map_type()) {
  if (!nice) {
    out << "[";
    for (size_t i = 0; i < list.size(); i++) {
      if (i > 0) {
        out << ", ";
      }
      fmt(out, list[i]);
    }
    out << "]";
  } else {
    out << "[\n";
    for (size_t i = 0; i < list.size(); i++) {
      if (i > 0) {
        out << ",\n";
      }
      out << "    ";
      fmt(out, list[i]);
    }
    out << "\n]\n";
  }
  return out;
}

template <typename Ctx>
static void print_hist(std::ostream& out, Ctx& ctx, size_t levels) {
  for (size_t level = 0; level < levels + 1; level++) {
    auto hist_size = ctx.hist_size(level);
    out << "hist[" << level << "]: [";
    for (size_t i = 0; i < hist_size; i++) {
      auto hist = ctx.hist(level, i);
      out << hist << ", ";
    }
    out << "]\n";
  }
}

/******************************************************************************/
