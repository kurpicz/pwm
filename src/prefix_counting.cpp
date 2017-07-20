#include "algorithm.hpp"
#include "wavelet_structures.hpp"
#include "prefix_counting.hpp"

void wt_pc_uint8_uint32(const uint8_t* text, const uint32_t size,
  const uint32_t levels) {
  construction::prefix_counting<uint8_t, uint32_t, wavelet_tree>(text, size, levels);
}

void wt_pc_uint16_uint32(const uint16_t* text, const uint32_t size,
  const uint32_t levels) {
  construction::prefix_counting<uint16_t, uint32_t, wavelet_tree>(text, size, levels);
}

void wt_pc_uint32_uint32(const uint32_t* text, const uint32_t size,
  const uint32_t levels) {
  construction::prefix_counting<uint32_t, uint32_t, wavelet_tree>(text, size, levels);
}

void wm_pc_uint8_uint32(const uint8_t* text, const uint32_t size,
  const uint32_t levels) {
  construction::prefix_counting<uint8_t, uint32_t, wavelet_matrix>(text, size, levels);
}

void wm_pc_uint16_uint32(const uint16_t* text, const uint32_t size,
  const uint32_t levels) {
  construction::prefix_counting<uint16_t, uint32_t, wavelet_matrix>(text, size, levels);
}

void wm_pc_uint32_uint32(const uint32_t* text, const uint32_t size,
  const uint32_t levels) {
  construction::prefix_counting<uint32_t, uint32_t, wavelet_matrix>(text, size, levels);
}

CONSTRUCTION_ALGORITHM_REGISTER("TEST", "TEST", uint8_t,uint32_t, wavelet_tree,
  wt_pc_uint8_uint32)