/*******************************************************************************
 * include/external_memory/internal_memory_bound.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

struct internal_memory_bound {
  static uint64_t& value() {
    static uint64_t bound = 4ULL * 1024 * 1024 * 1024;
    return bound;
  }
};

/******************************************************************************/