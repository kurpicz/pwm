/*******************************************************************************
 * include/wm_pddps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WM_DOMAIN_DECOMPOSITION_PREFIX_SORTING_PARALLEL
#define WM_DOMAIN_DECOMPOSITION_PREFIX_SORTING_PARALLEL

#include <cstring>
#include <omp.h>
#include <vector>

#include "common.hpp"

template <typename TextType, typename SizeType>
class wm_pddps {

public:
  wm_pddps(const std::vector<TextType>& text, const SizeType size,
    const SizeType levels) : _bv(levels), _zeros(levels, 0) {

  }

  auto get_bv_and_zeros() const {
    return std::make_pair(_bv, _zeros);
  }

private:
  std::vector<uint64_t*> _bv;
  std::vector<SizeType> _zeros;
}; // class wm_pddps

#endif // WM_DOMAIN_DECOMPOSITION_PREFIX_SORTING_PARALLEL

/******************************************************************************/