/*******************************************************************************
 * include/huffman/huff_m_level_sizes_builder.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/common.hpp"

template<typename AlphabetType>
class level_sizes_builder {
  histogram<AlphabetType> m_hist;
  std::vector<uint64_t> m_level_sizes;
public:
  // init
  level_sizes_builder(histogram<AlphabetType>&& h):
    m_hist(std::move(h)) {}
  inline void allocate_levels(size_t levels) {
    m_level_sizes = std::vector<uint64_t>(levels, 0);
  }
  inline void count(size_t level, AlphabetType symbol) {
    m_level_sizes[level] += m_hist.frequency(symbol);
  }
  inline histogram<AlphabetType> const& get_histogram() {
    return m_hist;
  }

  // finalization
  inline void drop_hist_and_finalize_level_sizes() {
    for (uint64_t i = m_level_sizes.size() - 1; i > 0; --i) {
      m_level_sizes[i - 1] += m_level_sizes[i];
    }
    drop_me(std::move(m_hist));
  }

  // final data
  inline size_t levels() const { return level_sizes().size(); }
  inline std::vector<uint64_t> const& level_sizes() const {
    return m_level_sizes;
  }
};

/******************************************************************************/
