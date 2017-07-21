#include "benchmark/algorithm.hpp"
#include "prefix_sorting.hpp"

using wm_ps_8 = wm_ps<uint8_t>;
using wm_ps_16 = wm_ps<uint16_t>;
using wm_ps_32 = wm_ps<uint32_t>;

CONSTRUCTION_REGISTER("wm_ps<uint8_t, uint32_t>",
  "Sqeuential wavelet matrix construction with 8-bit alphabet (using sorting).",
  wm_ps_8)
CONSTRUCTION_REGISTER("wm_ps<uint16_t, uint32_t>",
  "Sqeuential wavelet matrix construction with 16-bit alphabet (using sorting).",
  wm_ps_16)
CONSTRUCTION_REGISTER("wm_ps<uint32_t, uint32_t>",
  "Sqeuential wavelet matrix construction with 32-bit alphabet (using sorting).",
  wm_ps_32)