/*******************************************************************************
 * src/dd_prefix_counting.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "wx_dd_pc.hpp"
#include "benchmark/algorithm.hpp"

using wm_dd_pc_8 = wx_dd_pc<uint8_t, true>;
using wm_dd_pc_16 = wx_dd_pc<uint16_t, true>;
using wm_dd_pc_32 = wx_dd_pc<uint32_t, true>;

CONSTRUCTION_REGISTER("wm_dd_pc",
  "Parallel wavelet matrix construction with 8-bit alphabet "
  "(using domain decomposition and prefix counting).", wm_dd_pc_8)
CONSTRUCTION_REGISTER("wm_dd_pc",
  "Parallel wavelet matrix construction with 16-bit alphabet "
  "(using domain decomposition and prefix counting).", wm_dd_pc_16)
CONSTRUCTION_REGISTER("wm_dd_pc",
  "Parallel wavelet matrix construction with 32-bit alphabet "
  "(using domain decomposition and prefix counting).", wm_dd_pc_32)

using wt_dd_pc_8 = wx_dd_pc<uint8_t, false>;
using wt_dd_pc_16 = wx_dd_pc<uint16_t, false>;
using wt_dd_pc_32 = wx_dd_pc<uint32_t, false>;

CONSTRUCTION_REGISTER("wt_dd_pc",
  "Parallel wavelet tree construction with 8-bit alphabet "
  "(using domain decomposition and prefix counting).", wt_dd_pc_8)
CONSTRUCTION_REGISTER("wt_dd_pc",
  "Parallel wavelet tree construction with 16-bit alphabet "
  "(using domain decomposition and prefix counting).", wt_dd_pc_16)
CONSTRUCTION_REGISTER("wt_dd_pc",
  "Parallel wavelet tree construction with 32-bit alphabet "
  "(using domain decomposition and prefix counting).", wt_dd_pc_32)

/******************************************************************************/
