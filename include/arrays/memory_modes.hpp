/*******************************************************************************
 * include/util/memory_types.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

enum memory_mode : bool { internal = false, external = true };

//template <bool value>
//struct boolToMemMode {
//  static constexpr memory_mode result =
//      value ? memory_mode::external : memory_mode::internal;
//};

/******************************************************************************/
