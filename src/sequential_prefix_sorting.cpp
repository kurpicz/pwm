  /*******************************************************************************
 * src/sequential_prefix_sorting.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "wx_ps.hpp"
#include "wx_ps_oe.hpp"
#include "wx_ps_fe.hpp"
#include "benchmark/algorithm.hpp"

using wm_ps_8 = wx_ps<uint8_t, false, memory_mode::internal>;
using wm_ps_16 = wx_ps<uint16_t, false, memory_mode::internal>;
using wm_ps_32 = wx_ps<uint32_t, false, memory_mode::internal>;

using wt_ps_8 = wx_ps<uint8_t, true, memory_mode::internal>;
using wt_ps_16 = wx_ps<uint16_t, true, memory_mode::internal>;
using wt_ps_32 = wx_ps<uint32_t, true, memory_mode::internal>;

CONSTRUCTION_REGISTER("wm_ps",
  "Sequential wavelet matrix construction with 8-bit alphabet (using sorting).",
  wm_ps_8)
CONSTRUCTION_REGISTER("wm_ps",
  "Sequential wavelet matrix construction with 16-bit alphabet (using sorting).",
  wm_ps_16)
CONSTRUCTION_REGISTER("wm_ps",
  "Sequential wavelet matrix construction with 32-bit alphabet (using sorting).",
  wm_ps_32)

CONSTRUCTION_REGISTER("wt_ps",
  "Sequential wavelet tree construction with 8-bit alphabet (using sorting).",
  wt_ps_8)
CONSTRUCTION_REGISTER("wt_ps",
  "Sequential wavelet tree construction with 16-bit alphabet (using sorting).",
  wt_ps_16)
CONSTRUCTION_REGISTER("wt_ps",
  "Sequential wavelet tree construction with 32-bit alphabet (using sorting).",
  wt_ps_32)



using wm_ps_8_oe = wx_ps_oe<uint8_t, false, memory_mode::external_output>;
using wm_ps_16_oe = wx_ps_oe<uint16_t, false, memory_mode::external_output>;
using wm_ps_32_oe = wx_ps_oe<uint32_t, false, memory_mode::external_output>;

using wt_ps_8_oe = wx_ps_oe<uint8_t, true, memory_mode::external_output>;
using wt_ps_16_oe = wx_ps_oe<uint16_t, true, memory_mode::external_output>;
using wt_ps_32_oe = wx_ps_oe<uint32_t, true, memory_mode::external_output>;

CONSTRUCTION_REGISTER("wm_ps_oe",
  "Sequential wavelet matrix construction with 8-bit alphabet (using sorting, external output).",
  wm_ps_8_oe)
CONSTRUCTION_REGISTER("wm_ps_oe",
  "Sequential wavelet matrix construction with 16-bit alphabet (using sorting, external output).",
  wm_ps_16_oe)
CONSTRUCTION_REGISTER("wm_ps_oe",
  "Sequential wavelet matrix construction with 32-bit alphabet (using sorting, external output).",
  wm_ps_32_oe)

CONSTRUCTION_REGISTER("wt_ps_oe",
  "Sequential wavelet tree construction with 8-bit alphabet (using sorting, external output).",
  wt_ps_8_oe)
CONSTRUCTION_REGISTER("wt_ps_oe",
  "Sequential wavelet tree construction with 16-bit alphabet (using sorting, external output).",
  wt_ps_16_oe)
CONSTRUCTION_REGISTER("wt_ps_oe",
  "Sequential wavelet tree construction with 32-bit alphabet (using sorting, external output).",
  wt_ps_32_oe)



//~ using wm_ps_8_fe = wx_ps_fe<uint8_t, false>;
//~ using wm_ps_16_fe = wx_ps_fe<uint16_t, false>;
//~ using wm_ps_32_fe = wx_ps_fe<uint32_t, false>;

using wt_ps_8_fe = wx_ps_fe<uint8_t, true>;
using wt_ps_16_fe = wx_ps_fe<uint16_t, true>;
using wt_ps_32_fe = wx_ps_fe<uint32_t, true>;

//~ CONSTRUCTION_REGISTER("wm_ps_fe",
  //~ "Sequential wavelet matrix construction with 8-bit alphabet (using sorting, fully external).",
  //~ wm_ps_8_fe)
//~ CONSTRUCTION_REGISTER("wm_ps_fe",
  //~ "Sequential wavelet matrix construction with 16-bit alphabet (using sorting, fully external).",
  //~ wm_ps_16_fe)
//~ CONSTRUCTION_REGISTER("wm_ps_fe",
  //~ "Sequential wavelet matrix construction with 32-bit alphabet (using sorting, fully external).",
  //~ wm_ps_32_fe)

CONSTRUCTION_REGISTER("wt_ps_fe",
  "Sequential wavelet tree construction with 8-bit alphabet (using sorting, fully external).",
  wt_ps_8_fe)
CONSTRUCTION_REGISTER("wt_ps_fe",
  "Sequential wavelet tree construction with 16-bit alphabet (using sorting, fully external).",
  wt_ps_16_fe)
CONSTRUCTION_REGISTER("wt_ps_fe",
  "Sequential wavelet tree construction with 32-bit alphabet (using sorting, fully external).",
  wt_ps_32_fe)

/******************************************************************************/
