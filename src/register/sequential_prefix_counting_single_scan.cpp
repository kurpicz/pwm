/*******************************************************************************
 * src/register/sequential_prefix_counting.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wx_pc_ss.hpp"

using wm_pc_ss_8 = wx_pc_ss<uint8_t, false>;
using wm_pc_ss_16 = wx_pc_ss<uint16_t, false>;
using wm_pc_ss_32 = wx_pc_ss<uint32_t, false>;

CONSTRUCTION_REGISTER("wm_pc_ss",
  "Sequential wavelet matrix construction with 8-bit alphabet "
  "(using counting and computing all levels at once).", wm_pc_ss_8)
CONSTRUCTION_REGISTER("wm_pc_ss",
  "Sequential wavelet matrix construction with 16-bit alphabet "
  "(using counting and computing all levels at once).", wm_pc_ss_16)
CONSTRUCTION_REGISTER("wm_pc_ss",
  "Sequential wavelet matrix construction with 32-bit alphabet "
  "(using counting and computing all levels at once).", wm_pc_ss_32)

using wt_pc_ss_8 = wx_pc_ss<uint8_t, true>;
using wt_pc_ss_16 = wx_pc_ss<uint16_t, true>;
using wt_pc_ss_32 = wx_pc_ss<uint32_t, true>;

CONSTRUCTION_REGISTER("wt_pc_ss",
  "Sequential wavelet tree construction with 8-bit alphabet "
  "(using counting and computing all levels at once).", wt_pc_ss_8)
CONSTRUCTION_REGISTER("wt_pc_ss",
  "Sequential wavelet tree construction with 16-bit alphabet "
  "(using counting and computing all levels at once).", wt_pc_ss_16)
CONSTRUCTION_REGISTER("wt_pc_ss",
  "Sequential wavelet tree construction with 32-bit alphabet "
  "(using counting and computing all levels at once).", wt_pc_ss_32)

/******************************************************************************/
