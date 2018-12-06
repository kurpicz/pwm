/*******************************************************************************
 * src/register/sequential_prefix_counting.cpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wx_dd_fe.hpp"

template <uint64_t bytesMemory>
struct wx_dd_fe_registry {
  using wm_dd_fe = wx_dd_fe<uint8_t, false, bytesMemory>;
  using wt_dd_fe = wx_dd_fe<uint8_t, true, bytesMemory>;


  CONSTRUCTION_REGISTER_MEMBER(
      "wm_dd_fe",
      "Sequential wavelet matrix construction with 8-bit alphabet "
      "(using counting, external input).",
      wm_dd_fe)
  CONSTRUCTION_REGISTER_MEMBER(
      "wt_dd_fe",
      "Sequential wavelet tree construction with 8-bit alphabet "
      "(using counting, external input).",
      wt_dd_fe)
};

/******************************************************************************/
