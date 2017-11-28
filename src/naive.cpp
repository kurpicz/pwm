/*******************************************************************************
 * src/naive.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "huffman/wx_huff_naive.hpp"
#include "wx_naive.hpp"

using wm_naive_8 = wx_naive<uint8_t, false>;
using wm_naive_16 = wx_naive<uint16_t, false>;
using wm_naive_32 = wx_naive<uint32_t, false>;

CONSTRUCTION_REGISTER("wm_naive",
 "Naive sequential wavelet matrix construction with 8-bit alphabet.",
 wm_naive_8)
CONSTRUCTION_REGISTER("wm_naive",
  "Naive sequential wavelet matrix construction with 16-bit alphabet.",
  wm_naive_16)
CONSTRUCTION_REGISTER("wm_naive",
  "Naive sequential wavelet matrix construction with 32-bit alphabet.",
  wm_naive_32)

using wt_naive_8 = wx_naive<uint8_t, true>;
using wt_naive_16 = wx_naive<uint16_t, true>;
using wt_naive_32 = wx_naive<uint32_t, true>;

CONSTRUCTION_REGISTER("wt_naive",
  "Naive sequential wavelet tree construction with 8-bit alphabet.",
  wt_naive_8)
CONSTRUCTION_REGISTER("wt_naive",
  "Naive sequential wavelet tree construction with 16-bit alphabet.",
  wt_naive_16)
CONSTRUCTION_REGISTER("wt_naive",
  "Naive sequential wavelet tree construction with 32-bit alphabet.",
  wt_naive_32)

using wm_huff_naive_8 = wx_huff_naive<uint8_t, true>;
CONSTRUCTION_REGISTER("wm_huff_naive_8",
  "Naive sequential Huffman-shaped wavelet matrix construction with 8-bit"
  "alphabet.", wm_huff_naive_8)

using wt_huff_naive_8 = wx_huff_naive<uint8_t, false>;
CONSTRUCTION_REGISTER("wt_huff_naive_8",
  "Naive sequential Huffman-shaped wavelet tree construction with 8-bit"
  "alphabet.", wt_huff_naive_8)

/******************************************************************************/
