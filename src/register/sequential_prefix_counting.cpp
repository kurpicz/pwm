/*******************************************************************************
 * src/register/sequential_prefix_counting.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "util/type_for_bytes.hpp"
#include "benchmark/algorithm.hpp"
#include "wx_pc.hpp"

using wm_pc_8 = wx_pc<uint8_t, false>;
using wm_pc_16 = wx_pc<uint16_t, false>;
using wm_pc_32 = wx_pc<uint32_t, false>;
using wm_pc_40 = wx_pc<uint40_t, false>;

CONSTRUCTION_REGISTER("wm_pc",
  "Sequential wavelet matrix construction with 8-bit alphabet "
  "(using counting).", wm_pc_8)
CONSTRUCTION_REGISTER("wm_pc",
  "Sequential wavelet matrix construction with 16-bit alphabet "
  "(using counting).", wm_pc_16)
CONSTRUCTION_REGISTER("wm_pc",
  "Sequential wavelet matrix construction with 32-bit alphabet "
  "(using counting).", wm_pc_32)
CONSTRUCTION_REGISTER("wm_pc",
  "Sequential wavelet matrix construction with 40-bit alphabet "
  "(using counting).", wm_pc_40)

using wt_pc_8 = wx_pc<uint8_t, true>;
using wt_pc_16 = wx_pc<uint16_t, true>;
using wt_pc_32 = wx_pc<uint32_t, true>;
using wt_pc_40 = wx_pc<uint40_t, true>;

CONSTRUCTION_REGISTER("wt_pc",
  "Sequential wavelet tree construction with 8-bit alphabet "
  "(using counting).", wt_pc_8)
CONSTRUCTION_REGISTER("wt_pc",
  "Sequential wavelet tree construction with 16-bit alphabet "
  "(using counting).", wt_pc_16)
CONSTRUCTION_REGISTER("wt_pc",
  "Sequential wavelet tree construction with 32-bit alphabet "
  "(using counting).", wt_pc_32)
CONSTRUCTION_REGISTER("wt_pc",
  "Sequential wavelet tree construction with 40-bit alphabet "
  "(using counting).", wt_pc_40)

/******************************************************************************/
