/*******************************************************************************
 * src/parallel_prefix_sorting.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wx_pps.hpp"

using wm_pps_8 = wx_pps<uint8_t, false>;
using wm_pps_16 = wx_pps<uint16_t, false>;
using wm_pps_32 = wx_pps<uint32_t, false>;

CONSTRUCTION_REGISTER("wm_pps",
  "Parallel wavelet matrix construction with 8-bit alphabet "
  "(using sorting).", wm_pps_8, false, false)
CONSTRUCTION_REGISTER("wm_pps",
  "Parallel wavelet matrix construction with 16-bit alphabet "
  "(using sorting).", wm_pps_16, false, false)
CONSTRUCTION_REGISTER("wm_pps",
  "Parallel wavelet matrix construction with 32-bit alphabet "
  "(using sorting).", wm_pps_32, false, false)

using wt_pps_8 = wx_pps<uint8_t, true>;
using wt_pps_16 = wx_pps<uint16_t, true>;
using wt_pps_32 = wx_pps<uint32_t, true>;

CONSTRUCTION_REGISTER("wt_pps",
  "Parallel wavelet tree construction with 8-bit alphabet "
  "(using sorting).", wt_pps_8, false, false)
CONSTRUCTION_REGISTER("wt_pps",
  "Parallel wavelet tree construction with 16-bit alphabet "
  "(using sorting).", wt_pps_16, false, false)
CONSTRUCTION_REGISTER("wt_pps",
  "Parallel wavelet tree construction with 32-bit alphabet "
  "(using sorting).", wt_pps_32, false, false)

/******************************************************************************/
