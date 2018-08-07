/*******************************************************************************
 * src/parallel_prefix_counting.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wx_ppc.hpp"

using wm_ppc_8  = wx_ppc<uint8_t, false, memory_mode::internal>;
using wm_ppc_16 = wx_ppc<uint16_t, false, memory_mode::internal>;
using wm_ppc_32 = wx_ppc<uint32_t, false, memory_mode::internal>;

CONSTRUCTION_REGISTER("wm_ppc",
  "Parallel wavelet matrix construction with 8-bit alphabet "
  "(using counting).", wm_ppc_8)
CONSTRUCTION_REGISTER("wm_ppc",
  "Parallel wavelet matrix construction with 16-bit alphabet "
  "(using counting).", wm_ppc_16)
CONSTRUCTION_REGISTER("wm_ppc",
  "Parallel wavelet matrix construction with 32-bit alphabet "
  "(using counting).", wm_ppc_32)

using wt_ppc_8  = wx_ppc<uint8_t, true, memory_mode::internal>;
using wt_ppc_16 = wx_ppc<uint16_t, true, memory_mode::internal>;
using wt_ppc_32 = wx_ppc<uint32_t, true, memory_mode::internal>;

CONSTRUCTION_REGISTER("wt_ppc",
  "Parallel wavelet tree construction with 8-bit alphabet "
  "(using counting).", wt_ppc_8)
CONSTRUCTION_REGISTER("wt_ppc",
  "Parallel wavelet tree construction with 16-bit alphabet "
  "(using counting).", wt_ppc_16)
CONSTRUCTION_REGISTER("wt_ppc",
  "Parallel wavelet tree construction with 32-bit alphabet "
  "(using counting).", wt_ppc_32)

/******************************************************************************/
