  /*******************************************************************************
 * src/sequential_prefix_sorting.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "wx_ps.hpp"
#include "benchmark/algorithm.hpp"

using wm_ps_8 = wx_ps<uint8_t, false>;
using wm_ps_16 = wx_ps<uint16_t, false>;
using wm_ps_32 = wx_ps<uint32_t, false>;

CONSTRUCTION_REGISTER("wm_ps",
  "Sequential wavelet matrix construction with 8-bit alphabet (using sorting).",
  wm_ps_8, false, false)
CONSTRUCTION_REGISTER("wm_ps",
  "Sequential wavelet matrix construction with 16-bit alphabet (using sorting).",
  wm_ps_16, false, false)
CONSTRUCTION_REGISTER("wm_ps",
  "Sequential wavelet matrix construction with 32-bit alphabet (using sorting).",
  wm_ps_32, false, false)

using wt_ps_8 = wx_ps<uint8_t, true>;
using wt_ps_16 = wx_ps<uint16_t, true>;
using wt_ps_32 = wx_ps<uint32_t, true>;

CONSTRUCTION_REGISTER("wt_ps",
  "Sequential wavelet tree construction with 8-bit alphabet (using sorting).",
  wt_ps_8, false, false)
CONSTRUCTION_REGISTER("wt_ps",
  "Sequential wavelet tree construction with 16-bit alphabet (using sorting).",
  wt_ps_16, false, false)
CONSTRUCTION_REGISTER("wt_ps",
  "Sequential wavelet tree construction with 32-bit alphabet (using sorting).",
  wt_ps_32, false, false)


// external output

CONSTRUCTION_REGISTER("wm_ps_se",
  "Sequential wavelet matrix construction with 8-bit alphabet (using sorting, semi-external).",
  wm_ps_8, false, true)
CONSTRUCTION_REGISTER("wm_ps_se",
  "Sequential wavelet matrix construction with 16-bit alphabet (using sorting, semi-external).",
  wm_ps_16, false, true)
CONSTRUCTION_REGISTER("wm_ps_se",
  "Sequential wavelet matrix construction with 32-bit alphabet (using sorting, semi-external).",
  wm_ps_32, false, true)

CONSTRUCTION_REGISTER("wt_ps_se",
  "Sequential wavelet tree construction with 8-bit alphabet (using sorting, semi-external).",
  wt_ps_8, false, true)
CONSTRUCTION_REGISTER("wt_ps_se",
  "Sequential wavelet tree construction with 16-bit alphabet (using sorting, semi-external).",
  wt_ps_16, false, true)
CONSTRUCTION_REGISTER("wt_ps_se",
  "Sequential wavelet tree construction with 32-bit alphabet (using sorting, semi-external).",
  wt_ps_32, false, true)

/******************************************************************************/
