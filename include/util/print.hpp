/*******************************************************************************
 * include/util/debug.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <iostream>
#include <string>
#include <functional>

#include "construction/wavelet_structure.hpp"
#include "util/debug.hpp"

[[maybe_unused]]
static void print_structure(std::ostream& out, wavelet_structure const& structure) {
  const base_bit_vectors& bv = structure.bvs();
  const std::vector<uint64_t>& zeros = structure.zeros();
  for (uint64_t i = 0; i < bv.levels(); i++) {
    out << "   bv["<<i<<"]";
    out << "[";
    for (uint64_t j = 0; j < bv.level_bit_size(i); j++) {
      out << uint64_t(bit_at(bv[i], j)) << "";
    }
    out << "]";
    if (!structure.is_tree() && i != (bv.levels() - 1)) {
      out << " zeros[" << i << "] = " << zeros.at(i);
    }
    out << std::endl;
  }
  out << std::endl;
}

template<typename T>
struct force_integer_trait {
  inline static void print(std::ostream& out, T const& v) {
    out << v;
  }
};
template<>
struct force_integer_trait<uint8_t> {
  inline static void print(std::ostream& out, uint8_t const& v) {
    out << int(v);
}
};
/// Print T, but convert it to int first if it is a char-like type
struct print_force_type {
    template<typename T>
    inline void operator()(std::ostream& out, T const& t) const {
      force_integer_trait<T>::print(out, t);
    }
};

template<typename list_type, typename map_type = print_force_type>
static std::ostream& print_list(std::ostream& out, list_type const& list,
                                bool nice = false,
                                map_type fmt = map_type()) {
    if (!nice) {
        out << "[";
        for(size_t i = 0; i < list.size(); i++) {
            if (i > 0) {
                out << ", ";
            }
            fmt(out, list[i]);
        }
        out << "]";
    } else {
        out << "[\n";
        for(size_t i = 0; i < list.size(); i++) {
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

/******************************************************************************/
