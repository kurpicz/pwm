#pragma once

#include "common.hpp"

template<typename SizeType>
struct SingleThreaded {
    std::vector<SizeType> m_hist;
    std::vector<SizeType> m_bit_reverse;
    std::vector<SizeType> m_borders;

    SingleThreaded(SizeType const levels) {
        auto cur_max_char = (1ull << levels);

        m_hist.reserve(cur_max_char);
        m_hist.resize(cur_max_char, 0);

        m_bit_reverse = BitReverse<SizeType>(levels - 1);

        m_borders.reserve(cur_max_char);
        m_borders.resize(cur_max_char, 0);
    }

    SizeType const hist_size(SizeType const level) {
        return 1ull << level;
    }

    SizeType& hist(SizeType const level, SizeType const i) {
        return m_hist[i];
    }

    SizeType& rho(size_t level, size_t i) {
        return m_bit_reverse[i];
    }

    std::vector<SizeType>& borders() {
        return m_borders;
    }

    template<typename F>
    void if_zeros(F f) {
        f();
    }
};

template<typename SizeType>
struct DdMultiThreaded {
    DdMultiThreaded(SizeType const levels) {
    }

    SizeType const hist_size(SizeType const level) {
        return 1ull << level;
    }

    SizeType& hist(SizeType const level, SizeType const i) {

    }

    SizeType& rho(size_t level, size_t i) {

    }

    std::vector<SizeType>& borders() {

    }

    template<typename F>
    void if_zeros(F f) {

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
