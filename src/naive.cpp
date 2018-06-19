/*******************************************************************************
 * src/naive.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wx_naive.hpp"

using wm_naive_8 = wx_naive<uint8_t, false>;
using wm_naive_16 = wx_naive<uint16_t, false>;
using wm_naive_32 = wx_naive<uint32_t, false>;

CONSTRUCTION_REGISTER("wm_naive",
  "Sequential wavelet matrix construction with 8-bit alphabet "
  "(using counting).", wm_naive_8, false, false)
CONSTRUCTION_REGISTER("wm_naive",
  "Sequential wavelet matrix construction with 16-bit alphabet "
  "(using counting).", wm_naive_16, false, false)
CONSTRUCTION_REGISTER("wm_naive",
  "Sequential wavelet matrix construction with 32-bit alphabet "
  "(using counting).", wm_naive_32, false, false)

using wt_naive_8 = wx_naive<uint8_t, true>;
using wt_naive_16 = wx_naive<uint16_t, true>;
using wt_naive_32 = wx_naive<uint32_t, true>;

CONSTRUCTION_REGISTER("wt_naive",
  "Sequential wavelet tree construction with 8-bit alphabet "
  "(using counting).", wt_naive_8, false, false)
CONSTRUCTION_REGISTER("wt_naive",
  "Sequential wavelet tree construction with 16-bit alphabet "
  "(using counting).", wt_naive_16, false, false)
CONSTRUCTION_REGISTER("wt_naive",
  "Sequential wavelet tree construction with 32-bit alphabet "
  "(using counting).", wt_naive_32, false, false)

/******************************************************************************/
