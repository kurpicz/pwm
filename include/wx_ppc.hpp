#pragma once

#include "util/common.hpp"
#include "util/pc.hpp"
#include "util/wavelet_structure.hpp"

template <typename AlphabetType, bool is_matrix>
class wx_ppc {
  using ctx_t = KeepLevel<is_matrix,
    decltype(rho_dispatch<is_matrix>::create(0))>;

public:
  static constexpr bool    is_parallel = true;
  static constexpr bool    is_tree     = !is_matrix;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);

  static wavelet_structure compute(AlphabetType const* const global_text,
    const uint64_t size, const uint64_t levels) {

    if (size == 0) {
      return wavelet_structure();
    }

    const uint64_t omp_size = omp_get_max_threads();
  }
}; // class wx_ppc