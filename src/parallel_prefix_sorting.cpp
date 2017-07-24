/*******************************************************************************
 * src/sequential_prefix_sorting.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wm_pps.hpp"
#include "wt_pps.hpp"

using wm_pps_8 = wm_pps<uint8_t>;
using wm_pps_16 = wm_pps<uint16_t>;
using wm_pps_32 = wm_pps<uint32_t>;

CONSTRUCTION_REGISTER("wm_pps<uint8_t>",
  "Parallel wavelet matrix construction with 8-bit alphabet "
  "(using counting).", wm_pps_8)
CONSTRUCTION_REGISTER("wm_pps<uint16_t>",
  "Parallel wavelet matrix construction with 16-bit alphabet "
  "(using counting).", wm_pps_16)
CONSTRUCTION_REGISTER("wm_pps<uint32_t>",
  "Parallel wavelet matrix construction with 32-bit alphabet "
  "(using counting).", wm_pps_32)

using wt_pps_8 = wt_pps<uint8_t>;
using wt_pps_16 = wt_pps<uint16_t>;
using wt_pps_32 = wt_pps<uint32_t>;

CONSTRUCTION_REGISTER("wt_pps<uint8_t>",
  "Parallel wavelet tree construction with 8-bit alphabet "
  "(using counting).", wt_pps_8)
CONSTRUCTION_REGISTER("wt_pps<uint16_t>",
  "Parallel wavelet tree construction with 16-bit alphabet "
  "(using counting).", wt_pps_16)
CONSTRUCTION_REGISTER("wt_pps<uint32_t>",
  "Parallel wavelet tree construction with 32-bit alphabet "
  "(using counting).", wt_pps_32)

/******************************************************************************/
