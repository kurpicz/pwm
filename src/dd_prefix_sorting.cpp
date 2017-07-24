/*******************************************************************************
 * src/dd_prefix_sorting.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "wx_dd_ps.hpp"
#include "benchmark/algorithm.hpp"

using wm_dd_ps_8 = wx_dd_ps<uint8_t, true>;
using wm_dd_ps_16 = wx_dd_ps<uint16_t, true>;
using wm_dd_ps_32 = wx_dd_ps<uint32_t, true>;

CONSTRUCTION_REGISTER("wm_dd_ps<uint8_t>",
  "Parallel wavelet matrix construction with 8-bit alphabet "
  "(using domain decomposition and prefix sorting).", wm_dd_ps_8)
CONSTRUCTION_REGISTER("wm_dd_ps_<uint16_t>",
  "Parallel wavelet matrix construction with 16-bit alphabet "
  "(using domain decomposition and prefix sorting).", wm_dd_ps_16)
CONSTRUCTION_REGISTER("wm_dd_ps<uint32_t>",
  "Parallel wavelet matrix construction with 32-bit alphabet "
  "(using domain decomposition and prefix sorting).", wm_dd_ps_32)

using wt_dd_ps_8 = wx_dd_ps<uint8_t, false>;
using wt_dd_ps_16 = wx_dd_ps<uint16_t, false>;
using wt_dd_ps_32 = wx_dd_ps<uint32_t, false>;

CONSTRUCTION_REGISTER("wt_dd_ps<uint8_t>",
  "Parallel wavelet tree construction with 8-bit alphabet "
  "(using domain decomposition and prefix sorting).", wt_dd_ps_8)
CONSTRUCTION_REGISTER("wt_dd_ps_<uint16_t>",
  "Parallel wavelet tree construction with 16-bit alphabet "
  "(using domain decomposition and prefix sorting).", wt_dd_ps_16)
CONSTRUCTION_REGISTER("wt_dd_ps<uint32_t>",
  "Parallel wavelet tree construction with 32-bit alphabet "
  "(using domain decomposition and prefix sorting).", wt_dd_ps_32)

/******************************************************************************/
