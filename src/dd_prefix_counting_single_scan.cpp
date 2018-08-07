/*******************************************************************************
 * src/dd_prefix_counting_single_scan.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "wx_dd_pc_ss.hpp"
#include "benchmark/algorithm.hpp"

using wm_dd_pc_ss_8 = wx_dd_pc_ss<uint8_t, false, memory_mode::internal>;
using wm_dd_pc_ss_16 = wx_dd_pc_ss<uint16_t, false, memory_mode::internal>;
using wm_dd_pc_ss_32 = wx_dd_pc_ss<uint32_t, false, memory_mode::internal>;

CONSTRUCTION_REGISTER("wm_dd_pc_ss",
  "Parallel wavelet matrix construction with 8-bit alphabet "
  "(using domain decomposition and prefix counting and a single scan).",
  wm_dd_pc_ss_8)
CONSTRUCTION_REGISTER("wm_dd_pc_ss",
  "Parallel wavelet matrix construction with 16-bit alphabet "
  "(using domain decomposition and prefix counting and a single scan).",
  wm_dd_pc_ss_16)
CONSTRUCTION_REGISTER("wm_dd_pc_ss",
  "Parallel wavelet matrix construction with 32-bit alphabet "
  "(using domain decomposition and prefix counting and a single scan).",
  wm_dd_pc_ss_32)

using wt_dd_pc_ss_8 = wx_dd_pc_ss<uint8_t, true, memory_mode::internal>;
using wt_dd_pc_ss_16 = wx_dd_pc_ss<uint16_t, true, memory_mode::internal>;
using wt_dd_pc_ss_32 = wx_dd_pc_ss<uint32_t, true, memory_mode::internal>;

CONSTRUCTION_REGISTER("wt_dd_pc_ss",
  "Parallel wavelet tree construction with 8-bit alphabet "
  "(using domain decomposition and prefix counting and a single scan).",
  wt_dd_pc_ss_8)
CONSTRUCTION_REGISTER("wt_dd_pc_ss",
  "Parallel wavelet tree construction with 16-bit alphabet "
  "(using domain decomposition and prefix counting and a single scan).",
  wt_dd_pc_ss_16)
CONSTRUCTION_REGISTER("wt_dd_pc_ss",
  "Parallel wavelet tree construction with 32-bit alphabet "
  "(using domain decomposition and prefix counting and a single scan).",
  wt_dd_pc_ss_32)

/******************************************************************************/
