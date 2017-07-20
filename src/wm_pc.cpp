#include "algorithm.hpp"
#include "wm_pc.hpp"

using wm_pc_8 = wm_pc<uint8_t>;
using wm_pc_16 = wm_pc<uint16_t>;
using wm_pc_32 = wm_pc<uint32_t>;
CONSTRUCTION_REGISTER("wm_pc<uint8_t, uint32_t>",
  "Sqeuential wavelet matrix construction with 8-bit alphabet.", wm_pc_8)
CONSTRUCTION_REGISTER("wm_pc<uint16_t, uint32_t>",
  "Sqeuential wavelet matrix construction with 16-bit alphabet.", wm_pc_16)
CONSTRUCTION_REGISTER("wm_pc<uint32_t, uint32_t>",
  "Sqeuential wavelet matrix construction with 32-bit alphabet.", wm_pc_32)