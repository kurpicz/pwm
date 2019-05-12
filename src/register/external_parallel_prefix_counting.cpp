/*******************************************************************************
 * src/register/sequential_prefix_counting.cpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "external_memory/wx_ppc_ie.hpp"

using wm_ppc_8_ie = wx_ppc_ie<uint8_t, false>;
using wt_ppc_8_ie = wx_ppc_ie<uint8_t, true>;

CONSTRUCTION_REGISTER(
    "wm_ppc_ie",
    "Parallel wavelet matrix construction with 8-bit alphabet "
    "(using counting, external input).",
    wm_ppc_8_ie)
CONSTRUCTION_REGISTER(
    "wt_ppc_ie",
    "Parallel wavelet tree construction with 8-bit alphabet "
    "(using counting, external input).",
    wt_ppc_8_ie)

/******************************************************************************/
