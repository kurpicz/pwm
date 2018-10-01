/*******************************************************************************
 * src/register/huff_naive.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wx_huff_naive.hpp"

using wm_huff_naive_8 = wx_huff_naive<uint8_t, false>;
using wm_huff_naive_16 = wx_huff_naive<uint16_t, false>;
using wm_huff_naive_32 = wx_huff_naive<uint32_t, false>;
CONSTRUCTION_REGISTER("wm_huff_naive",
  "Naive sequential Huffman-shaped wavelet matrix construction with 8-bit"
  " alphabet.", wm_huff_naive_8)
CONSTRUCTION_REGISTER("wm_huff_naive",
  "Naive sequential Huffman-shaped wavelet matrix construction with 16-bit"
  " alphabet.", wm_huff_naive_16)
CONSTRUCTION_REGISTER("wm_huff_naive",
  "Naive sequential Huffman-shaped wavelet matrix construction with 32-bit"
  " alphabet.", wm_huff_naive_32)

using wt_huff_naive_8 = wx_huff_naive<uint8_t, true>;
using wt_huff_naive_16 = wx_huff_naive<uint16_t, true>;
using wt_huff_naive_32 = wx_huff_naive<uint32_t, true>;
CONSTRUCTION_REGISTER("wt_huff_naive",
  "Naive sequential Huffman-shaped wavelet tree construction with 8-bit"
  " alphabet.", wt_huff_naive_8)
CONSTRUCTION_REGISTER("wt_huff_naive",
  "Naive sequential Huffman-shaped wavelet tree construction with 16-bit"
  " alphabet.", wt_huff_naive_16)
CONSTRUCTION_REGISTER("wt_huff_naive",
  "Naive sequential Huffman-shaped wavelet tree construction with 32-bit"
  " alphabet.", wt_huff_naive_32)

/******************************************************************************/
