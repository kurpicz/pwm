/*******************************************************************************
 * src/register/sequential_prefix_counting.cpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "external_memory/wx_dd_pc_fe.hpp"

using wt_dd_pc_fe1 = wx_dd_pc_fe<uint8_t, true, false>;
using wt_dd_pc_fe2 = wx_dd_pc_fe<uint8_t, true, true>;

using wm_dd_pc_fe1 = wx_dd_pc_fe<uint8_t, false, false>;
using wm_dd_pc_fe2 = wx_dd_pc_fe<uint8_t, false, true>;

CONSTRUCTION_REGISTER_MEMBER(
    "wt_dd_pc_fe1",
    "Parallel wavelet tree construction with 8-bit alphabet "
    "(using domain decomposition and prefix counting, fully external).",
    wt_dd_pc_fe1)

CONSTRUCTION_REGISTER_MEMBER(
    "wt_dd_pc_fe2",
    "Parallel wavelet tree construction with 8-bit alphabet "
    "(using domain decomposition and prefix counting, fully external).",
    wt_dd_pc_fe2)

CONSTRUCTION_REGISTER_MEMBER(
    "wm_dd_pc_fe1",
    "Parallel wavelet matrix construction with 8-bit alphabet "
    "(using domain decomposition and prefix counting, fully external).",
    wm_dd_pc_fe1)

CONSTRUCTION_REGISTER_MEMBER(
    "wm_dd_pc_fe2",
    "Parallel wavelet matrix construction with 8-bit alphabet "
    "(using domain decomposition and prefix counting, fully external).",
    wm_dd_pc_fe2)

/******************************************************************************/
