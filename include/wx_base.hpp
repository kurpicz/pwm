/*******************************************************************************
 * include/wx_base.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/memory_types.hpp"

template <
  typename AlphabetType, 
  bool is_tree_, 
  bool is_parallel_, 
  bool is_huffman_shaped_, 
  memory_mode mem_mode_>
class wx_base {

public:
  static constexpr bool  is_parallel = is_parallel_;
  static constexpr bool  is_tree   = is_tree_;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);
  static constexpr bool  is_huffman_shaped = is_huffman_shaped_;
  static constexpr memory_mode mem_mode = mem_mode_;
}; // class wx_ps


#define WX_BASE(AlphabetType, is_tree_, is_parallel_, is_huffman_shaped_, mem_mode_) \
  static constexpr bool  is_parallel = is_parallel_; \
  static constexpr bool  is_tree   = is_tree_; \
  static constexpr uint8_t word_width  = sizeof(AlphabetType); \
  static constexpr bool  is_huffman_shaped = is_huffman_shaped_; \
  static constexpr memory_mode mem_mode = mem_mode_; \

/******************************************************************************/
