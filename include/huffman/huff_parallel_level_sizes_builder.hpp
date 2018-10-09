/*******************************************************************************
 * include/huffman/huff_parallel_level_sizes_builder.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/common.hpp"

template <typename AlphabetType>
class parallel_level_sizes_builder {
  histogram<AlphabetType> m_hist;
  std::vector<histogram<AlphabetType>> m_hists;

  std::vector<uint64_t> global_level_sizes;
  std::vector<std::vector<uint64_t>> local_level_sizes;

  inline size_t shards() {
    return m_hists.size();
  }

public:
  // init
  parallel_level_sizes_builder(histogram<AlphabetType>&& h,
                               std::vector<histogram<AlphabetType>>&& hs)
      : m_hist(h), m_hists(hs) {}
  inline void allocate_levels(size_t levels) {
    global_level_sizes = std::vector<uint64_t>(levels, 0);
    local_level_sizes.resize(shards());
    for (auto& e : local_level_sizes) {
      e = std::vector<uint64_t>(levels, 0);
    }
  }
  inline void count(size_t level, AlphabetType symbol) {
    for (size_t i = 0; i < shards(); i++) {
      auto freq = m_hists[i].frequency(symbol);
      global_level_sizes[level] += freq;
      local_level_sizes[i][level] += freq;
    }
  }
  inline histogram<AlphabetType> const& get_histogram() const {
    return m_hist;
  }

  // finalization
  inline void drop_hist_and_finalize_level_sizes() {
    for (uint64_t i = global_level_sizes.size() - 1; i > 0; --i) {
      global_level_sizes[i - 1] += global_level_sizes[i];
    }
    for (size_t j = 0; j < shards(); j++) {
      for (uint64_t i = local_level_sizes[j].size() - 1; i > 0; --i) {
        local_level_sizes[j][i - 1] += local_level_sizes[j][i];
      }
    }
    drop_me(std::move(m_hist));
    drop_me(std::move(m_hists));
  }

  // final data
  inline size_t levels() const {
    return level_sizes().size();
  }
  inline std::vector<uint64_t> const& level_sizes() const {
    return global_level_sizes;
  }
  inline std::vector<std::vector<uint64_t>> const& level_sizes_shards() const {
    return local_level_sizes;
  }
};

/******************************************************************************/
