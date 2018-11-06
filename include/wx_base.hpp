/*******************************************************************************
 * include/wx_base.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "arrays/memory_modes.hpp"



//TODO: Add other attributes

template <bool external_in_, bool external_out_>
class wx_in_out_external {
public:
  static constexpr bool external_in = external_in_;
  static constexpr bool external_out = external_out_;
}; // class wx_ps



/******************************************************************************/
