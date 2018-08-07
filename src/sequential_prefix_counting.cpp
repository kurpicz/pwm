/*******************************************************************************
 * src/sequential_prefix_counting.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wx_pc.hpp"
#include "wx_pc_ie.hpp"


using wm_pc_8 = wx_pc<uint8_t, false, memory_mode::internal>;
using wm_pc_16 = wx_pc<uint16_t, false, memory_mode::internal>;
using wm_pc_32 = wx_pc<uint32_t, false, memory_mode::internal>;

using wt_pc_8 = wx_pc<uint8_t, true, memory_mode::internal>;
using wt_pc_16 = wx_pc<uint16_t, true, memory_mode::internal>;
using wt_pc_32 = wx_pc<uint32_t, true, memory_mode::internal>;

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

using wm_pc_8_ie = wx_pc_ie<uint8_t, false, memory_mode::external_input>;
using wm_pc_16_ie = wx_pc_ie<uint16_t, false, memory_mode::external_input>;
using wm_pc_32_ie = wx_pc_ie<uint32_t, false, memory_mode::external_input>;

using wt_pc_8_ie = wx_pc_ie<uint8_t, true, memory_mode::external_input>;
using wt_pc_16_ie = wx_pc_ie<uint16_t, true, memory_mode::external_input>;
using wt_pc_32_ie = wx_pc_ie<uint32_t, true, memory_mode::external_input>;


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


//~ using wm_pc_8_fe = wx_pc<uint8_t, false, memory_mode::external>;
//~ using wm_pc_16_fe = wx_pc<uint16_t, false, memory_mode::external>;
//~ using wm_pc_32_fe = wx_pc<uint32_t, false, memory_mode::external>;

//~ using wt_pc_8_fe = wx_pc<uint8_t, true, memory_mode::external>;
//~ using wt_pc_16_fe = wx_pc<uint16_t, true, memory_mode::external>;
//~ using wt_pc_32_fe = wx_pc<uint32_t, true, memory_mode::external>;

//~ CONSTRUCTION_REGISTER("wm_pc_fe",
  //~ "Sequential wavelet matrix construction with 8-bit alphabet "
  //~ "(using counting, external).", wm_pc_8_fe)
//~ CONSTRUCTION_REGISTER("wm_pc_fe",
  //~ "Sequential wavelet matrix construction with 16-bit alphabet "
  //~ "(using counting, external).", wm_pc_16_fe)
//~ CONSTRUCTION_REGISTER("wm_pc_fe",
  //~ "Sequential wavelet matrix construction with 32-bit alphabet "
  //~ "(using counting, external).", wm_pc_32_fe)

//~ CONSTRUCTION_REGISTER("wt_pc_fe",
  //~ "Sequential wavelet tree construction with 8-bit alphabet "
  //~ "(using counting, external).", wt_pc_8_fe)
//~ CONSTRUCTION_REGISTER("wt_pc_fe",
  //~ "Sequential wavelet tree construction with 16-bit alphabet "
  //~ "(using counting, external).", wt_pc_16_fe)
//~ CONSTRUCTION_REGISTER("wt_pc_fe",
  //~ "Sequential wavelet tree construction with 32-bit alphabet "
  //~ "(using counting, external).", wt_pc_32_fe)

/******************************************************************************/
