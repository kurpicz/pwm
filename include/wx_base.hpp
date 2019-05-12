/*******************************************************************************
 * include/wx_base.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

// TODO: Add other attributes

template <bool external_in_,
          bool external_out_,
          bool stats_support_ = false,
          bool is_inplace_ = false>
class wx_in_out_external {
public:
  static constexpr bool external_in = external_in_;
  static constexpr bool external_out = external_out_;
  static constexpr bool stats_support = stats_support_;
  static constexpr bool is_inplace = is_inplace_;
}; // class wx_ps

/******************************************************************************/
