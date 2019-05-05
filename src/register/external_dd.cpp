/*******************************************************************************
 * src/register/sequential_prefix_counting.cpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "wx_dd_pc_fe.hpp"
#include "construction/pc_dd_fe/wx_dd_pc_fe2.hpp"

template <uint64_t bytesMemory>
struct wx_dd_fe_registry {
//  using wm_dd_fe = wx_dd_fe<uint8_t, false, bytesMemory>;
//  using wt_dd_fe = wx_dd_fe<uint8_t, true, bytesMemory>;

//  using wm_dd_pc_fe_0 = wx_dd_pc_fe2<uint8_t, false, bytesMemory>;
  using wt_dd_pc_feNEW_1 = wx_dd_pc_fe2<uint8_t, true, bytesMemory, false>;
  using wt_dd_pc_feNEW_2 = wx_dd_pc_fe2<uint8_t, true, bytesMemory, true>;

//  using wm_dd_pc_fe_1 = wx_dd_pc_fe<uint8_t, false, bytesMemory, false>;
  using wt_dd_pc_fe_1 = wx_dd_pc_fe<uint8_t, true, bytesMemory / 8, false>;
  using wt_dd_pc_fe_2 = wx_dd_pc_fe<uint8_t, true, bytesMemory / 8, true>;


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

//  CONSTRUCTION_REGISTER_MEMBER(
//      "wm_dd_pc_fe_0",
//      "Parallel wavelet matrix construction with 8-bit alphabet "
//      "(using domain decomposition and prefix counting, fully external).",
//      wm_dd_pc_fe_0)
//  CONSTRUCTION_REGISTER_MEMBER(
//      "wt_dd_pc_fe_NEW",
//      "Parallel wavelet tree construction with 8-bit alphabet "
//      "(using domain decomposition and prefix counting, fully external).",
//      wt_dd_pc_fe_0)

//  CONSTRUCTION_REGISTER_MEMBER(
//      "wm_dd_pc_fe_1",
//      "Parallel wavelet matrix construction with 8-bit alphabet "
//      "(using domain decomposition and prefix counting, fully external).",
//      wm_dd_pc_fe_1)

//  CONSTRUCTION_REGISTER_MEMBER(
//      "wt_dd_pc_feNEW_1",
//      "Parallel wavelet tree construction with 8-bit alphabet "
//      "(using domain decomposition and prefix counting, fully external).",
//      wt_dd_pc_feNEW_1)

  CONSTRUCTION_REGISTER_MEMBER(
      "wt_dd_pc_feNEW_2",
      "Parallel wavelet tree construction with 8-bit alphabet "
      "(using domain decomposition and prefix counting, fully external).",
      wt_dd_pc_feNEW_2)

  CONSTRUCTION_REGISTER_MEMBER(
      "wt_dd_pc_fe_1",
      "Parallel wavelet tree construction with 8-bit alphabet "
      "(using domain decomposition and prefix counting, fully external).",
      wt_dd_pc_fe_1)

//  CONSTRUCTION_REGISTER_MEMBER(
//      "wt_dd_pc_fe_2",
//      "Parallel wavelet tree construction with 8-bit alphabet "
//      "(using domain decomposition and prefix counting, fully external).",
//      wt_dd_pc_fe_2)
};

/******************************************************************************/
