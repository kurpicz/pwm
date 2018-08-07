/*******************************************************************************
 * src/parallel_prefix_counting_single_scan.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wx_ppc_ss.hpp"

using wm_ppc_ss_8 = wx_ppc_ss<uint8_t, false, memory_mode::internal>;
using wm_ppc_ss_16 = wx_ppc_ss<uint16_t, false, memory_mode::internal>;
using wm_ppc_ss_32 = wx_ppc_ss<uint32_t, false, memory_mode::internal>;

CONSTRUCTION_REGISTER("wm_ppc_ss",
  "Parallel wavelet matrix construction with 8-bit alphabet "
  "(using counting and computing all levels at once).", wm_ppc_ss_8)
CONSTRUCTION_REGISTER("wm_ppc_ss",
  "Parallel wavelet matrix construction with 16-bit alphabet "
  "(using counting and computing all levels at once).", wm_ppc_ss_16)
CONSTRUCTION_REGISTER("wm_ppc_ss",
  "Parallel wavelet matrix construction with 32-bit alphabet "
  "(using counting and computing all levels at once).", wm_ppc_ss_32)

using wt_ppc_ss_8 = wx_ppc_ss<uint8_t, true, memory_mode::internal>;
using wt_ppc_ss_16 = wx_ppc_ss<uint16_t, true, memory_mode::internal>;
using wt_ppc_ss_32 = wx_ppc_ss<uint32_t, true, memory_mode::internal>;

CONSTRUCTION_REGISTER("wt_ppc_ss",
  "Parallel wavelet tree construction with 8-bit alphabet "
  "(using counting and computing all levels at once).", wt_ppc_ss_8)
CONSTRUCTION_REGISTER("wt_ppc_ss",
  "Parallel wavelet tree construction with 16-bit alphabet "
  "(using counting and computing all levels at once).", wt_ppc_ss_16)
CONSTRUCTION_REGISTER("wt_ppc_ss",
  "Parallel wavelet tree construction with 32-bit alphabet "
  "(using counting and computing all levels at once).", wt_ppc_ss_32)

/******************************************************************************/
