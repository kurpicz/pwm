#pragma once

#include "common.hpp"

template<typename SizeType>
struct SingleThreaded {
    SizeType m_cur_max_char;
    std::vector<SizeType> m_hist;

    SingleThreaded(SizeType const levels, SizeType const cur_max_char):
        m_cur_max_char(cur_max_char)
    {
        m_hist.reserve(cur_max_char);
        m_hist.resize(cur_max_char, 0);
    }

    SizeType const hist_size(SizeType const level) {
        return 1ull << level;
    }

    SizeType& hist(SizeType const level, SizeType const i) {
        return m_hist[i];
    }

    SizeType& cur_max_char() {
        return m_cur_max_char;
    }
};

template<typename SizeType>
struct DdMultiThreaded {
    SizeType m_cur_max_char;

    DdMultiThreaded(SizeType const levels, SizeType const cur_max_char) {
        m_cur_max_char = cur_max_char;
    }

    SizeType const hist_size(SizeType const level) {
        return 1ull << level;
    }

    SizeType& hist(SizeType const level, SizeType const i) {

    }

    SizeType& cur_max_char() {
        return m_cur_max_char;
    }
};


template<
    typename AlphabetType,
    typename SizeType,
    typename Context
>
Bvs<SizeType> pc(AlphabetType const* text,
                 SizeType const size,
                 SizeType const levels,
                 Context& context)
{
    auto bv = Bvs<SizeType>(size, levels);



    return bv;
}
