#include "benchmark/algorithm.hpp"
#include "wx_pc.hpp"

using wm_pc_8 = wx_pc<uint8_t, true>;
using wm_pc_16 = wx_pc<uint16_t, true>;
using wm_pc_32 = wx_pc<uint32_t, true>;

CONSTRUCTION_REGISTER("wm_pc<uint8_t, uint32_t>",
  "Sqeuential wavelet matrix construction with 8-bit alphabet (using counting).", wm_pc_8)
CONSTRUCTION_REGISTER("wm_pc<uint16_t, uint32_t>",
  "Sqeuential wavelet matrix construction with 16-bit alphabet (using counting).", wm_pc_16)
CONSTRUCTION_REGISTER("wm_pc<uint32_t, uint32_t>",
  "Sqeuential wavelet matrix construction with 32-bit alphabet (using counting).", wm_pc_32)