/*******************************************************************************
 * src/register/huff_ps.cpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wx_huff_ps.hpp"

using wm_huff_ps_8 = wx_huff_ps<uint8_t, false>;
using wm_huff_ps_16 = wx_huff_ps<uint16_t, false>;
using wm_huff_ps_32 = wx_huff_ps<uint32_t, false>;
CONSTRUCTION_REGISTER("huff_wm_ps",
  "Sequential Huffman-shaped wavelet matrix construction with 8-bit"
  " alphabet (using sorting).", wm_huff_ps_8)
CONSTRUCTION_REGISTER("huff_wm_ps",
  "Sequential Huffman-shaped wavelet matrix construction with 16-bit"
  " alphabet (using sorting).", wm_huff_ps_16)
CONSTRUCTION_REGISTER("huff_wm_ps",
  "Sequential Huffman-shaped wavelet matrix construction with 32-bit"
  " alphabet (using sorting).", wm_huff_ps_32)

using wt_huff_ps_8 = wx_huff_ps<uint8_t, true>;
using wt_huff_ps_16 = wx_huff_ps<uint16_t, true>;
using wt_huff_ps_32 = wx_huff_ps<uint32_t, true>;
CONSTRUCTION_REGISTER("huff_wt_ps",
  "Sequential Huffman-shaped wavelet tree construction with 8-bit"
  " alphabet (using sorting).", wt_huff_ps_8)
CONSTRUCTION_REGISTER("huff_wt_ps",
  "Sequential Huffman-shaped wavelet tree construction with 16-bit"
  " alphabet (using sorting).", wt_huff_ps_16)
CONSTRUCTION_REGISTER("huff_wt_ps",
  "Sequential Huffman-shaped wavelet tree construction with 32-bit"
  " alphabet (using sorting).", wt_huff_ps_32)

/******************************************************************************/
