/*******************************************************************************
 * src/register/sequential_prefix_counting.cpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
//#include "wx_dd_fe.hpp"
#include "wx_dd_pc_fe.hpp"

template <uint64_t bytesMemory>
struct wx_dd_fe_registry {
//  using wm_dd_fe = wx_dd_fe<uint8_t, false, bytesMemory>;
//  using wt_dd_fe = wx_dd_fe<uint8_t, true, bytesMemory>;

  using wm_dd_pc_fe = wx_dd_pc_fe<uint8_t, false, bytesMemory>;
  using wt_dd_pc_fe = wx_dd_pc_fe<uint8_t, true, bytesMemory>;


//  CONSTRUCTION_REGISTER_MEMBER(
//      "wm_dd_fe",
//      "Parallel wavelet matrix construction with 8-bit alphabet "
//      "(using domain decomposition, fully external).",
//      wm_dd_fe)
//  CONSTRUCTION_REGISTER_MEMBER(
//      "wt_dd_fe",
//      "Parallel wavelet tree construction with 8-bit alphabet "
//      "(using domain decomposition, fully external).",
//      wt_dd_fe)

  CONSTRUCTION_REGISTER_MEMBER(
      "wm_dd_pc_fe",
      "Parallel wavelet matrix construction with 8-bit alphabet "
      "(using domain decomposition and prefix counting, fully external).",
      wm_dd_pc_fe)
  CONSTRUCTION_REGISTER_MEMBER(
      "wt_dd_pc_fe",
      "Parallel wavelet tree construction with 8-bit alphabet "
      "(using domain decomposition and prefix counting, fully external).",
      wt_dd_pc_fe)
};

/******************************************************************************/
