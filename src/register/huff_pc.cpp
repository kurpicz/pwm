/*******************************************************************************
 * src/register/huff_pc.cpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wx_huff_pc.hpp"

using wm_huff_pc_8 = wx_huff_pc<uint8_t, false>;
using wm_huff_pc_16 = wx_huff_pc<uint16_t, false>;
using wm_huff_pc_32 = wx_huff_pc<uint32_t, false>;
CONSTRUCTION_REGISTER("wm_huff_pc",
  "Sequential Huffman-shaped wavelet matrix construction with 8-bit"
  " alphabet (using counting).", wm_huff_pc_8)
CONSTRUCTION_REGISTER("wm_huff_pc",
  "Sequential Huffman-shaped wavelet matrix construction with 16-bit"
  " alphabet (using counting).", wm_huff_pc_16)
CONSTRUCTION_REGISTER("wm_huff_pc",
  "Sequential Huffman-shaped wavelet matrix construction with 32-bit"
  " alphabet (using counting).", wm_huff_pc_32)

using wt_huff_pc_8 = wx_huff_pc<uint8_t, true>;
using wt_huff_pc_16 = wx_huff_pc<uint16_t, true>;
using wt_huff_pc_32 = wx_huff_pc<uint32_t, true>;
CONSTRUCTION_REGISTER("wt_huff_pc",
  "Sequential Huffman-shaped wavelet tree construction with 8-bit"
  " alphabet (using counting).", wt_huff_pc_8)
CONSTRUCTION_REGISTER("wt_huff_pc",
  "Sequential Huffman-shaped wavelet tree construction with 16-bit"
  " alphabet (using counting).", wt_huff_pc_16)
CONSTRUCTION_REGISTER("wt_huff_pc",
  "Sequential Huffman-shaped wavelet tree construction with 32-bit"
  " alphabet (using counting).", wt_huff_pc_32)

/******************************************************************************/
