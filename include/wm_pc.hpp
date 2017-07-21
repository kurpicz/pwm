/*******************************************************************************
 * include/wm_pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WM_PREFIX_COUNTING
#define WM_PREFIX_COUNTING

#include <vector>

#include "common.hpp"
#include "pc.hpp"

template <typename AlphabetType, typename SizeType = uint64_t>
class wm_pc {

public:
    wm_pc() = default;

    wm_pc(const std::vector<AlphabetType>& text,
          const SizeType size,
          const SizeType levels) {

        if(text.size() == 0) { return; }

        auto ctx = LevelSinglePass<SizeType, true>(size, levels);

        pc(text, size, levels, ctx);

        _zeros = std::move(ctx.zeros());
        _bv = std::move(ctx.bv());
    }

    auto get_bv_and_zeros() const {
        return std::make_pair(_bv.vec(), _zeros);
    }

private:
    Bvs<SizeType> _bv;
    std::vector<SizeType> _zeros;
}; // class wm_pc

#endif // WM_PREFIX_COUNTING

/******************************************************************************/
