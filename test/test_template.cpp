#include <gtest/gtest.h>
#include "test/util.hpp"

#include <unordered_map>

// WX_NAME is defined by the cmake file to be one of the algorithms
#include "@WX_NAME@.hpp"
#define WX_NAME @WX_NAME@
#cmakedefine WT_TYPE
#cmakedefine WM_TYPE

#include "common.hpp"

template<typename bv_t>
auto bit_at(const bv_t& bv, size_t i) -> uint8_t {
    size_t offset = i / 64ull;
    size_t word_offset = i % 64ull;
    return (bv[offset] >> (63ull - word_offset)) & 1ull;
}

template <typename AlphabetType>
auto ConstructWX(std::vector<AlphabetType>& text, const bool already_reduced) {
    uint64_t max_char = 0;
    if (!already_reduced) {
        std::unordered_map<AlphabetType, uint64_t> word_list;
        for (const AlphabetType& character : text) {
        auto result = word_list.find(character);
        if (result == word_list.end()) {
            word_list.emplace(character, max_char++);
        }
        }
        --max_char;
        for (uint64_t i = 0; i < text.size(); ++i) {
        text[i] = static_cast<AlphabetType>(word_list.find(text[i])->second);
        }
    } else {
        for (const AlphabetType& character : text) {
        if (character > max_char) {
            max_char = character;
        }
        }
    }
    //std::cout << "Greatest character: " << static_cast<uint64_t>(max_char) << std::endl;

    AlphabetType levels = 0;
    while (max_char) {
        max_char >>= 1;
        ++levels;
    }
    //std::cout << "Levels in WM: " << static_cast<uint64_t>(levels) << std::endl;

    return WX_NAME<AlphabetType, uint32_t>(text, text.size(), levels);
}

std::vector<std::vector<size_t>> level_sizes(
    const std::vector<uint64_t*> bv,
    uint64_t bit_offset,
    uint64_t bit_length,
    size_t level
) {
    if (level == bv.size()) {
        return {};
    }

    uint64_t zeroes = 0;

    for(uint64_t i = 0; i < bit_length; i++) {
        uint8_t bit = bit_at(bv[level], bit_offset + i);

        if (bit == 0) {
            zeroes++;
        }
    }

    uint64_t size_left = zeroes;
    uint64_t size_right = bit_length - zeroes;

    auto sizes_left  = level_sizes(bv, bit_offset,              size_left, level + 1);
    auto sizes_right = level_sizes(bv, bit_offset + size_left, size_right, level + 1);

    std::vector<std::vector<size_t>> r;

    r.push_back({});
    r.back().push_back({bit_length});
    //r.back().push_back({size_right});

    for (size_t j = 0; j < sizes_left.size(); j++) {
        r.push_back({});
        auto& v = r.back();

        for(auto& e : sizes_left[j]) {
            v.push_back(e);
        }
        for(auto& e : sizes_right[j]) {
            v.push_back(e);
        }
    }

    return r;
}

std::string decode_wt(const std::vector<uint64_t*> bv, size_t length) {
    auto ls = level_sizes(bv, 0, length, 0);

    for (auto& v : ls) {
        for (size_t i = 1; i < v.size(); i++) {
            v[i] = v[i - 1] + v[i];
        }
        for (size_t i = 1; i < v.size(); i++) {
            size_t j = v.size() - i;
            v[j] = v[j - 1];
        }
        if (v.size() > 0) {
            v[0] = 0;
        }
    }

    auto r = std::vector<uint8_t>(length);

    for (size_t i = 0; i < length; i++) {
        uint8_t value = 0;
        size_t j = 0;
        for (size_t level = 0; level < bv.size(); level++) {
            auto& offset = ls[level][j];
            uint8_t bit = bit_at(bv[level], offset);

            value <<= 1;
            value |= bit;

            offset++;
            j = 2 * j + bit;
        }
        r[i] = value;
    }

    return std::string(r.begin(), r.end());
}

std::string decode_wm(const std::vector<uint64_t*> bv,
                       const std::vector<uint32_t>& zeros,
                       size_t length) {
    if (bv.size() == 0) {
        return {};
    }

    /*
    for (size_t i = 0; i < bv.size(); i++) {
        std::cout << "   bv["<<i<<"]";

        std::cout << "[";
        for (size_t j = 0; j < length; j++) {
            std::cout << size_t(bit_at(bv[i], j)) << "";
        }
        std::cout << "], zeros: " << zeros[i];

        std::cout << "\n";
    }
    std::cout << "\n";
    */

    auto r = std::vector<uint8_t>(length);
    auto rtmp = std::vector<uint8_t>(length);

    for(size_t level = bv.size() - 1; level > 0; level--) {
        size_t offset0 = 0;
        size_t offset1 = zeros[level - 1];

        for(size_t i = 0; i < length; i++) {
            r[i] |= (bit_at(bv[level], i) << (bv.size() - level - 1));
        }

        for(size_t i = 0; i < length; i++) {
            if(bit_at(bv[level - 1], i) == 0) {
                rtmp[i] = r[offset0];
                offset0++;
            } else {
                rtmp[i] = r[offset1];
                offset1++;
            }
        }

        r.swap(rtmp);
    }

    for(size_t i = 0; i < length; i++) {
        r[i] |= bit_at(bv[0], i) << (bv.size() - 1);
    }

    return std::string(r.begin(), r.end());
}

TEST(WX_NAME, smoketest) {
    test::roundtrip_batch([](const std::string& s){
        auto vec = std::vector<uint8_t>(s.begin(), s.end());
        //std::cout << "string: '" << s << "'\n";
        //ConstructWX(vec, false);
        auto r = ConstructWX(vec, true);
#ifdef WT_TYPE
        auto bv = r.get_bv();
        auto decoded_s = decode_wt(bv, vec.size());
        ASSERT_EQ(s, decoded_s);
#endif
#ifdef WM_TYPE
        auto bvz = r.get_bv_and_zeros();
        auto decoded_s = decode_wm(bvz.first, bvz.second, vec.size());
        ASSERT_EQ(s, decoded_s);
#endif
    });
}
