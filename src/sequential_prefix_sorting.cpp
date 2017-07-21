#include "wx_ps.hpp"
#include "benchmark/algorithm.hpp"

using wm_ps_8 = wx_ps<uint8_t, true>;
using wm_ps_16 = wx_ps<uint16_t, true>;
using wm_ps_32 = wx_ps<uint32_t, true>;

CONSTRUCTION_REGISTER("wm_ps<uint8_t>",
  "Sqeuential wavelet matrix construction with 8-bit alphabet (using sorting).",
  wm_ps_8)
CONSTRUCTION_REGISTER("wm_ps<uint16_t>",
  "Sqeuential wavelet matrix construction with 16-bit alphabet (using sorting).",
  wm_ps_16)
CONSTRUCTION_REGISTER("wm_ps<uint32_t>",
  "Sqeuential wavelet matrix construction with 32-bit alphabet (using sorting).",
  wm_ps_32)

using wt_ps_8 = wx_ps<uint8_t, false>;
using wt_ps_16 = wx_ps<uint16_t, false>;
using wt_ps_32 = wx_ps<uint32_t, false>;

CONSTRUCTION_REGISTER("wm_ps<uint8_t>",
  "Sqeuential wavelet tree construction with 8-bit alphabet (using sorting).",
  wt_ps_8)
CONSTRUCTION_REGISTER("wm_ps<uint16_t>",
  "Sqeuential wavelet tree construction with 16-bit alphabet (using sorting).",
  wt_ps_16)
CONSTRUCTION_REGISTER("wm_ps<uint32_t>",
  "Sqeuential wavelet tree construction with 32-bit alphabet (using sorting).",
  wt_ps_32)