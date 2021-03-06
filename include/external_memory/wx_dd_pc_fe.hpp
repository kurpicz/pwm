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
#include "external_memory/full/dd/external_dd.hpp"
#include "construction/wavelet_structure.hpp"
#include "wavelet_structure_external.hpp"
#include "internal_memory_bound.hpp"


#include "wx_base.hpp"
#include "external_memory/wx_dd_pc_fe.hpp"

template <typename AlphabetType, bool is_tree_, bool rw_simultaneously>
class wx_dd_pc_fe : public wx_in_out_external<true, true, true> {
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

    WX_DD_PC_FE_VERBOSE << "Starting wx_dd_pc_fe. MiB memory: "
                        << (internal_memory_bound::value() / 1024ULL / 1024ULL)
                        << ".\n";

    stats.phase("dd");
    WX_DD_PC_FE_VERBOSE << "Initializing result...\n";
    // create empty result
    std::ostringstream name;
    name << "w" << (is_tree_ ? "t" : "m") << "_dd_fe";
    auto result =
        wavelet_structure_external_factory(is_tree_).
            histograms(is_tree_).zeros(!is_tree_).
            construct(size, levels, name.str(), 0);
    if(size == 0) return result;

    WX_DD_PC_FE_VERBOSE << "Creating context...\n";
    external_dd_ctx<InputType, rw_simultaneously> ctx(
        text, size, levels, internal_memory_bound::value());

    WX_DD_PC_FE_VERBOSE << "Running dd phase...\n";

    ctx.dd();

    WX_DD_PC_FE_VERBOSE << "Running merge phase...\n";
    stats.phase("merge");
    ctx.template merge<is_tree>(result);

    return result;
  }
}; // class wx_ps

/******************************************************************************/
