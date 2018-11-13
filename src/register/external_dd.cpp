/*******************************************************************************
 * src/register/sequential_prefix_counting.cpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wx_dd_fe.hpp"

using wm_dd_fe = wx_dd_fe<uint8_t, false>;
using wt_dd_fe = wx_dd_fe<uint8_t, true>;


CONSTRUCTION_REGISTER(
    "wm_dd_fe",
    "Sequential wavelet matrix construction with 8-bit alphabet "
    "(using counting, external input).",
    wm_dd_fe)
CONSTRUCTION_REGISTER(
    "wt_dd_fe",
    "Sequential wavelet tree construction with 8-bit alphabet "
    "(using counting, external input).",
    wt_dd_fe)


/******************************************************************************/
