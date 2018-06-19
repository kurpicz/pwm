/*******************************************************************************
 * src/sequential_prefix_counting.cpp
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
  "(using counting and computing all levels at once).", wm_pc_ss_8, false, false)
CONSTRUCTION_REGISTER("wm_pc_ss",
  "Sequential wavelet matrix construction with 16-bit alphabet "
  "(using counting and computing all levels at once).", wm_pc_ss_16, false, false)
CONSTRUCTION_REGISTER("wm_pc_ss",
  "Sequential wavelet matrix construction with 32-bit alphabet "
  "(using counting and computing all levels at once).", wm_pc_ss_32, false, false)

using wt_pc_ss_8 = wx_pc_ss<uint8_t, true>;
using wt_pc_ss_16 = wx_pc_ss<uint16_t, true>;
using wt_pc_ss_32 = wx_pc_ss<uint32_t, true>;

CONSTRUCTION_REGISTER("wt_pc_ss",
  "Sequential wavelet tree construction with 8-bit alphabet "
  "(using counting and computing all levels at once).", wt_pc_ss_8, false, false)
CONSTRUCTION_REGISTER("wt_pc_ss",
  "Sequential wavelet tree construction with 16-bit alphabet "
  "(using counting and computing all levels at once).", wt_pc_ss_16, false, false)
CONSTRUCTION_REGISTER("wt_pc_ss",
  "Sequential wavelet tree construction with 32-bit alphabet "
  "(using counting and computing all levels at once).", wt_pc_ss_32, false, false)

//~ using wm_pc_ss_8_se = wx_pc_ss<uint8_t, false, true>;
//~ using wm_pc_ss_16_se = wx_pc_ss<uint16_t, false, true>;
//~ using wm_pc_ss_32_se = wx_pc_ss<uint32_t, false, true>;

//~ CONSTRUCTION_REGISTER_SE("wm_pc_ss_se",
  //~ "Sequential wavelet matrix construction with 8-bit alphabet "
  //~ "(using counting and computing all levels at once, semi-external).", wm_pc_ss_8_se)
//~ CONSTRUCTION_REGISTER_SE("wm_pc_ss_se",
  //~ "Sequential wavelet matrix construction with 16-bit alphabet "
  //~ "(using counting and computing all levels at once, semi-external).", wm_pc_ss_16_se)
//~ CONSTRUCTION_REGISTER_SE("wm_pc_ss_se",
  //~ "Sequential wavelet matrix construction with 32-bit alphabet "
  //~ "(using counting and computing all levels at once, semi-external).", wm_pc_ss_32_se)

//~ using wt_pc_ss_8_se = wx_pc_ss<uint8_t, true, true>;
//~ using wt_pc_ss_16_se = wx_pc_ss<uint16_t, true, true>;
//~ using wt_pc_ss_32_se = wx_pc_ss<uint32_t, true, true>;

//~ CONSTRUCTION_REGISTER_SE("wt_pc_ss_se",
  //~ "Sequential wavelet tree construction with 8-bit alphabet "
  //~ "(using counting and computing all levels at once, semi-external).", wt_pc_ss_8_se)
//~ CONSTRUCTION_REGISTER_SE("wt_pc_ss_se",
  //~ "Sequential wavelet tree construction with 16-bit alphabet "
  //~ "(using counting and computing all levels at once, semi-external).", wt_pc_ss_16_se)
//~ CONSTRUCTION_REGISTER_SE("wt_pc_ss_se",
  //~ "Sequential wavelet tree construction with 32-bit alphabet "
  //~ "(using counting and computing all levels at once, semi-external).", wt_pc_ss_32_se)

/******************************************************************************/
