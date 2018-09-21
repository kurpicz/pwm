/*******************************************************************************
 * src/register/dd_prefix_sorting.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "wx_huff_dd_naive.hpp"
#include "benchmark/algorithm.hpp"

using wm_huff_dd_naive_8 = wx_huff_dd_naive<uint8_t, false>;
using wm_huff_dd_naive_16 = wx_huff_dd_naive<uint16_t, false>;
using wm_huff_dd_naive_32 = wx_huff_dd_naive<uint32_t, false>;

CONSTRUCTION_REGISTER("wm_huff_dd_naive",
  "Parallel wavelet matrix construction with 8-bit alphabet "
  "(using domain decomposition and the naive algorithm).", wm_huff_dd_naive_8)
CONSTRUCTION_REGISTER("wm_huff_dd_naive",
  "Parallel wavelet matrix construction with 16-bit alphabet "
  "(using domain decomposition and the naive algorithm).", wm_huff_dd_naive_16)
CONSTRUCTION_REGISTER("wm_huff_dd_naive",
  "Parallel wavelet matrix construction with 32-bit alphabet "
  "(using domain decomposition and the naive algorithm).", wm_huff_dd_naive_32)

using wt_huff_dd_naive_8 = wx_huff_dd_naive<uint8_t, true>;
using wt_huff_dd_naive_16 = wx_huff_dd_naive<uint16_t, true>;
using wt_huff_dd_naive_32 = wx_huff_dd_naive<uint32_t, true>;

CONSTRUCTION_REGISTER("wt_huff_dd_naive",
  "Parallel wavelet tree construction with 8-bit alphabet "
  "(using domain decomposition and the naive algorithm).", wt_huff_dd_naive_8)
CONSTRUCTION_REGISTER("wt_huff_dd_naive",
  "Parallel wavelet tree construction with 16-bit alphabet "
  "(using domain decomposition and the naive algorithm).", wt_huff_dd_naive_16)
CONSTRUCTION_REGISTER("wt_huff_dd_naive",
  "Parallel wavelet tree construction with 32-bit alphabet "
  "(using domain decomposition and the naive algorithm).", wt_huff_dd_naive_32)

/******************************************************************************/
