/*******************************************************************************
 * include/construction/pc_dd_fe/wx_dd_pc_fe2.hpp
 *
 * Copyright (C) 2019 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#define WX_DD_PC_FE_VERBOSE if constexpr (false) atomic_out()

#include <bitset>
#include <vector>
#include <thread>
#include <sstream>
#include <omp.h>

#include "util/print.hpp"
#include "arrays/memory_types.hpp"
#include "construction/pc_dd_fe/external_dd.hpp"
#include "construction/wavelet_structure.hpp"
#include "construction/wavelet_structure_external.hpp"


#include "wx_base.hpp"
#include "wx_dd_pc_fe.hpp"

template <typename AlphabetType, bool is_tree_, uint64_t bytes_memory_, bool rw_simultaneously>
class wx_dd_pc_fe2 : public wx_in_out_external<true, true, true> {
public:
  static constexpr bool is_parallel = true;
  static constexpr bool is_tree = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = false;

  template <typename InputType, typename stats_type>
  static wavelet_structure_external
  compute(const InputType& text,
          const uint64_t size,
          const uint64_t levels,
          stats_type& stats) {
    static_assert(sizeof(typename InputType::value_type) == sizeof(AlphabetType));

    stats.phase("dd");

    // create empty result
    std::ostringstream name;
    name << "w" << (is_tree_ ? "t" : "m") << "_dd_fe";
    auto result =
        wavelet_structure_external_factory(is_tree_).
            histograms(is_tree_).zeros(!is_tree_).
            construct(size, levels, name.str(), 0);
    if(size == 0) return result;

    external_dd_ctx<InputType, bytes_memory_, rw_simultaneously> ctx(text, size, levels);

    ctx.dd();
    stats.phase("merge");
    ctx.merge(result);

//    print_structure(std::cout, result.getInternalStructure(), true);
//
//    auto result2 = wx_dd_pc_fe<typename InputType::value_type, true, bytes_memory_ / 8>::compute(text, size, levels, stats);
//    print_structure(std::cout, result2.getInternalStructure(), true);

    return result;
  }
}; // class wx_ps

/******************************************************************************/
