/*******************************************************************************
 * include/wx_ps_fe.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>
#include "wx_base.hpp"
#include "util/ctx_single_level.hpp"
#include "util/wavelet_structure.hpp"
#include "util/ps.hpp"
#include "util/memory_types.hpp"

template <typename AlphabetType, bool is_tree_>
class wx_ps_fe {
public:

  WX_BASE(AlphabetType, is_tree_, false, false, memory_mode::external)

  template <typename InputType>
  static wavelet_structure compute(const InputType& text, const uint64_t size,
    const uint64_t levels) {

    using ctx_t = ctx_single_level<OutputType, wx_base<AlphabetType, is_tree_, false, false, mem_mode>::is_tree>;

    if(size == 0) { return wavelet_structure(); }

    auto ctx = ctx_t(size, levels);
    ps_fully_external<AlphabetType>(text, size, levels, ctx);

    //~ if (ctx_t::compute_zeros)  {
      //~ return wavelet_structure(std::move(ctx.bv()), std::move(ctx.zeros()));
    //~ } else {
      //~ return wavelet_structure(std::move(ctx.bv()));
    //~ }
    
    return wavelet_structure();
  }
}; // class wx_ps

/******************************************************************************/
