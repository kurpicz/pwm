/*******************************************************************************
 * src/register/sequential_prefix_counting.cpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
//#include "wx_ps_fe.hpp"
#include "wx_ps_oe.hpp"

// semi-external

using wm_ps_8_oe = wx_ps_oe<uint8_t, false>;
//using wm_ps_16_oe = wx_ps_oe<uint16_t, false>;
//using wm_ps_32_oe = wx_ps_oe<uint32_t, false>;
//
using wt_ps_8_oe = wx_ps_oe<uint8_t, true>;
//using wt_ps_16_oe = wx_ps_oe<uint16_t, true>;
//using wt_ps_32_oe = wx_ps_oe<uint32_t, true>;
//
CONSTRUCTION_REGISTER("wm_ps_oe",
                      "Sequential wavelet matrix construction with 8-bit "
                      "alphabet (using sorting, external output).",
                      wm_ps_8_oe)
//CONSTRUCTION_REGISTER("wm_ps_oe",
//                      "Sequential wavelet matrix construction with 16-bit "
//                      "alphabet (using sorting, external output).",
//                      wm_ps_16_oe)
//CONSTRUCTION_REGISTER("wm_ps_oe",
//                      "Sequential wavelet matrix construction with 32-bit "
//                      "alphabet (using sorting, external output).",
//                      wm_ps_32_oe)
//
CONSTRUCTION_REGISTER("wt_ps_oe",
                      "Sequential wavelet tree construction with 8-bit "
                      "alphabet (using sorting, external output).",
                      wt_ps_8_oe)
//CONSTRUCTION_REGISTER("wt_ps_oe",
//                      "Sequential wavelet tree construction with 16-bit "
//                      "alphabet (using sorting, external output).",
//                      wt_ps_16_oe)
//CONSTRUCTION_REGISTER("wt_ps_oe",
//                      "Sequential wavelet tree construction with 32-bit "
//                      "alphabet (using sorting, external output).",
//                      wt_ps_32_oe)


using wm_ps_8_oe_ip = wx_ps_oe<uint8_t, false, true>;
//using wm_ps_16_oe_ip = wx_ps_oe<uint16_t, false, true>;
//using wm_ps_32_oe_ip = wx_ps_oe<uint32_t, false, true>;
//
using wt_ps_8_oe_ip = wx_ps_oe<uint8_t, true, true>;
//using wt_ps_16_oe_ip = wx_ps_oe<uint16_t, true, true>;
//using wt_ps_32_oe_ip = wx_ps_oe<uint32_t, true, true>;

CONSTRUCTION_REGISTER("wm_ps_oe_inplace",
                      "Sequential wavelet matrix construction with 8-bit "
                      "alphabet (using inplace sorting, external output).",
                      wm_ps_8_oe_ip)
//CONSTRUCTION_REGISTER("wm_ps_oe_inplace",
//                      "Sequential wavelet matrix construction with 16-bit "
//                      "alphabet (using inplace sorting, external output).",
//                      wm_ps_16_oe_ip)
//CONSTRUCTION_REGISTER("wm_ps_oe_inplace",
//                      "Sequential wavelet matrix construction with 32-bit "
//                      "alphabet (using inplace sorting, external output).",
//                      wm_ps_32_oe_ip)
//
CONSTRUCTION_REGISTER("wt_ps_oe_inplace",
                      "Sequential wavelet tree construction with 8-bit "
                      "alphabet (using inplace sorting, external output).",
                      wt_ps_8_oe_ip)
//CONSTRUCTION_REGISTER("wt_ps_oe_inplace",
//                      "Sequential wavelet tree construction with 16-bit "
//                      "alphabet (using inplace sorting, external output).",
//                      wt_ps_16_oe_ip)
//CONSTRUCTION_REGISTER("wt_ps_oe_inplace",
//                      "Sequential wavelet tree construction with 32-bit "
//                      "alphabet (using inplace sorting, external output).",
//                      wt_ps_32_oe_ip)

// external

//using wm_ps_8_fe_wp0 = wx_ps_fe<uint8_t, false, 0>;
//using wt_ps_8_fe_wp0 = wx_ps_fe<uint8_t, true, 0>;
//
//using wm_ps_8_fe_wp1 = wx_ps_fe<uint8_t, false, 1>;
//using wt_ps_8_fe_wp1 = wx_ps_fe<uint8_t, true, 1>;
//
//using wm_ps_8_fe_wp2 = wx_ps_fe<uint8_t, false, 2>;
//using wt_ps_8_fe_wp2 = wx_ps_fe<uint8_t, true, 2>;
//
//CONSTRUCTION_REGISTER(
//    "wm_ps_fe_wp0",
//    "Sequential wavelet matrix construction with 8-bit alphabet (using "
//    "sorting, fully external, no wordpacking).",
//    wm_ps_8_fe_wp0)
//CONSTRUCTION_REGISTER(
//    "wt_ps_fe_wp0",
//    "Sequential wavelet tree construction with 8-bit alphabet (using sorting, "
//    "fully external, no wordpacking).",
//    wt_ps_8_fe_wp0)
//CONSTRUCTION_REGISTER(
//    "wm_ps_fe_wp1",
//    "Sequential wavelet matrix construction with 8-bit alphabet (using "
//    "sorting, fully external, wordpacking with padding).",
//    wm_ps_8_fe_wp1)
//CONSTRUCTION_REGISTER(
//    "wt_ps_fe_wp1",
//    "Sequential wavelet tree construction with 8-bit alphabet (using sorting, "
//    "fully external, wordpacking with padding).",
//    wt_ps_8_fe_wp1)
//CONSTRUCTION_REGISTER(
//    "wm_ps_fe_wp2",
//    "Sequential wavelet matrix construction with 8-bit alphabet (using "
//    "sorting, fully external, wordpacking without padding).",
//    wm_ps_8_fe_wp2)
//CONSTRUCTION_REGISTER(
//    "wt_ps_fe_wp2",
//    "Sequential wavelet tree construction with 8-bit alphabet (using sorting, "
//    "fully external, wordpacking without padding).",
//    wt_ps_8_fe_wp2)

/******************************************************************************/
