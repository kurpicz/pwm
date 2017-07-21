/*******************************************************************************
 * include/wt_prefix_counting.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WT_PREFIX_COUNTING
#define WT_PREFIX_COUNTING

#include <vector>

#include "common.hpp"
#include "pc.hpp"

template <typename AlphabetType, typename SizeType>
class wt_pc {

public:
    wt_pc() = delete;

    wt_pc(const std::vector<AlphabetType>& text,
          const SizeType size,
          const SizeType levels) {

        if(text.size() == 0) { return; }

        auto ctx = LevelSinglePass<SizeType, false>(size, levels);

        pc(text, size, levels, ctx);

        _bv = std::move(ctx.bv());
    }

    auto get_bv() const {
        return _bv.vec();
    }

private:
    Bvs<SizeType> _bv;
}; // class wt_pc

#endif // WT_PREFIX_COUNTING

/******************************************************************************/
