/*******************************************************************************
 * include/util/common.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <unistd.h>

inline static uint64_t get_cache_line_size() {
  int64_t cl = sysconf (_SC_LEVEL1_DCACHE_LINESIZE);
  return (cl > 0) ? (uint64_t)cl : 64;
}

//*****************************************************************************/
