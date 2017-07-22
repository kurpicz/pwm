/*******************************************************************************
 * src/sequential_prefix_counting.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wm_ppc.hpp"
#include "wt_ppc.hpp"

using wm_ppc_8 = wm_ppc<uint8_t>;
using wm_ppc_16 = wm_ppc<uint16_t>;
using wm_ppc_32 = wm_ppc<uint32_t>;

CONSTRUCTION_REGISTER("wm_ppc<uint8_t>",
  "Sqeuential wavelet matrix construction with 8-bit alphabet "
  "(using counting).", wm_ppc_8)
CONSTRUCTION_REGISTER("wm_ppc<uint16_t>",
  "Sqeuential wavelet matrix construction with 16-bit alphabet "
  "(using counting).", wm_ppc_16)
CONSTRUCTION_REGISTER("wm_ppc<uint32_t>",
  "Sqeuential wavelet matrix construction with 32-bit alphabet "
  "(using counting).", wm_ppc_32)

using wt_ppc_8 = wt_ppc<uint8_t>;
using wt_ppc_16 = wt_ppc<uint16_t>;
using wt_ppc_32 = wt_ppc<uint32_t>;

CONSTRUCTION_REGISTER("wt_ppc<uint8_t>",
  "Sqeuential wavelet matrix construction with 8-bit alphabet "
  "(using counting).", wt_ppc_8)
CONSTRUCTION_REGISTER("wt_ppc<uint16_t>",
  "Sqeuential wavelet matrix construction with 16-bit alphabet "
  "(using counting).", wt_ppc_16)
CONSTRUCTION_REGISTER("wt_ppc<uint32_t>",
  "Sqeuential wavelet matrix construction with 32-bit alphabet "
  "(using counting).", wt_ppc_32)

/******************************************************************************/
