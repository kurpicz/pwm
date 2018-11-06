/*******************************************************************************
 * src/register/huff_dd_pc_ss.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "wx_huff_dd_pc_ss.hpp"
#include "benchmark/algorithm.hpp"

using wm_huff_dd_pc_ss_8 = wx_huff_dd_pc_ss<uint8_t, false>;
using wm_huff_dd_pc_ss_16 = wx_huff_dd_pc_ss<uint16_t, false>;
using wm_huff_dd_pc_ss_32 = wx_huff_dd_pc_ss<uint32_t, false>;

CONSTRUCTION_REGISTER("wm_huff_dd_pc_ss",
  "Parallel wavelet matrix construction with 8-bit alphabet "
  "(using domain decomposition and the pc_ss algorithm).", wm_huff_dd_pc_ss_8)
CONSTRUCTION_REGISTER("wm_huff_dd_pc_ss",
  "Parallel wavelet matrix construction with 16-bit alphabet "
  "(using domain decomposition and the pc_ss algorithm).", wm_huff_dd_pc_ss_16)
CONSTRUCTION_REGISTER("wm_huff_dd_pc_ss",
  "Parallel wavelet matrix construction with 32-bit alphabet "
  "(using domain decomposition and the pc_ss algorithm).", wm_huff_dd_pc_ss_32)

using wt_huff_dd_pc_ss_8 = wx_huff_dd_pc_ss<uint8_t, true>;
using wt_huff_dd_pc_ss_16 = wx_huff_dd_pc_ss<uint16_t, true>;
using wt_huff_dd_pc_ss_32 = wx_huff_dd_pc_ss<uint32_t, true>;

CONSTRUCTION_REGISTER("wt_huff_dd_pc_ss",
  "Parallel wavelet tree construction with 8-bit alphabet "
  "(using domain decomposition and the pc_ss algorithm).", wt_huff_dd_pc_ss_8)
CONSTRUCTION_REGISTER("wt_huff_dd_pc_ss",
  "Parallel wavelet tree construction with 16-bit alphabet "
  "(using domain decomposition and the pc_ss algorithm).", wt_huff_dd_pc_ss_16)
CONSTRUCTION_REGISTER("wt_huff_dd_pc_ss",
  "Parallel wavelet tree construction with 32-bit alphabet "
  "(using domain decomposition and the pc_ss algorithm).", wt_huff_dd_pc_ss_32)

/******************************************************************************/
