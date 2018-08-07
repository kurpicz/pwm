/*******************************************************************************
 * include/wx_ps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
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

template <typename AlphabetType, bool is_tree_, memory_mode mem_mode_>
class wx_ps_oe {
public:

  WX_BASE(AlphabetType, is_tree_, false, false, mem_mode_)

  template <typename InputType, typename OutputType>
  static wavelet_structure<OutputType> compute(const InputType& text, const uint64_t size,
    const uint64_t levels) {

    using ctx_t = ctx_single_level<OutputType, wx_base<AlphabetType, is_tree_, false, false, mem_mode_>::is_tree>;

    if(size == 0) { return wavelet_structure<OutputType>(); }

    auto ctx = ctx_t(size, levels);

    auto sorted_text = std::vector<AlphabetType>(size);
    ps_out_external(text, size, levels, ctx, sorted_text.data());

    //~ if (ctx_t::compute_zeros)  {
      //~ return wavelet_structure<OutputType>(std::move(ctx.bv()), std::move(ctx.zeros()));
    //~ } else {
      //~ return wavelet_structure<OutputType>(std::move(ctx.bv()));
    //~ }
    
    return wavelet_structure<OutputType>();
  }
}; // class wx_ps

/******************************************************************************/
