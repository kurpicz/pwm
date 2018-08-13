/*******************************************************************************
 * src/sequential_prefix_counting.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wx_pc.hpp"
#include "wx_pc_ie.hpp"

using wm_pc_8 = wx_pc<uint8_t, false>;
using wm_pc_16 = wx_pc<uint16_t, false>;
using wm_pc_32 = wx_pc<uint32_t, false>;

using wt_pc_8 = wx_pc<uint8_t, true>;
using wt_pc_16 = wx_pc<uint16_t, true>;
using wt_pc_32 = wx_pc<uint32_t, true>;

CONSTRUCTION_REGISTER("wm_pc",
  "Sequential wavelet matrix construction with 8-bit alphabet "
  "(using counting).", wm_pc_8)
CONSTRUCTION_REGISTER("wm_pc",
  "Sequential wavelet matrix construction with 16-bit alphabet "
  "(using counting).", wm_pc_16)
CONSTRUCTION_REGISTER("wm_pc",
  "Sequential wavelet matrix construction with 32-bit alphabet "
  "(using counting).", wm_pc_32)

CONSTRUCTION_REGISTER("wt_pc",
  "Sequential wavelet tree construction with 8-bit alphabet "
  "(using counting).", wt_pc_8)
CONSTRUCTION_REGISTER("wt_pc",
  "Sequential wavelet tree construction with 16-bit alphabet "
  "(using counting).", wt_pc_16)
CONSTRUCTION_REGISTER("wt_pc",
  "Sequential wavelet tree construction with 32-bit alphabet "
  "(using counting).", wt_pc_32)

// semi-external

using wm_pc_8_ie = wx_pc_ie<uint8_t, false>;
using wm_pc_16_ie = wx_pc_ie<uint16_t, false>;
using wm_pc_32_ie = wx_pc_ie<uint32_t, false>;

using wt_pc_8_ie = wx_pc_ie<uint8_t, true>;
using wt_pc_16_ie = wx_pc_ie<uint16_t, true>;
using wt_pc_32_ie = wx_pc_ie<uint32_t, true>;

CONSTRUCTION_REGISTER("wm_pc_ie",
  "Sequential wavelet matrix construction with 8-bit alphabet "
  "(using counting, external input).", wm_pc_8_ie)
CONSTRUCTION_REGISTER("wm_pc_ie",
  "Sequential wavelet matrix construction with 16-bit alphabet "
  "(using counting, external input).", wm_pc_16_ie)
CONSTRUCTION_REGISTER("wm_pc_ie",
  "Sequential wavelet matrix construction with 32-bit alphabet "
  "(using counting, external input).", wm_pc_32_ie)

CONSTRUCTION_REGISTER("wt_pc_ie",
  "Sequential wavelet tree construction with 8-bit alphabet "
  "(using counting, external input).", wt_pc_8_ie)
CONSTRUCTION_REGISTER("wt_pc_ie",
  "Sequential wavelet tree construction with 16-bit alphabet "
  "(using counting, external input).", wt_pc_16_ie)
CONSTRUCTION_REGISTER("wt_pc_ie",
  "Sequential wavelet tree construction with 32-bit alphabet "
  "(using counting, external input).", wt_pc_32_ie)

/******************************************************************************/
