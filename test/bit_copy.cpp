#include <gtest/gtest.h>
#include "test/util.hpp"

#include "merge.hpp"

void bit_compare(
    std::vector<uint8_t> const& left,
    std::vector<uint8_t> const& right
) {
    auto left_s = bit_string<uint8_t>(left, left.size());
    auto right_s = bit_string<uint8_t>(right, right.size());

    if (left_s != right_s) {
        auto diff = std::string(std::max(left_s.size(), right_s.size()), ' ');

        size_t i = 0;
        for(; i < std::min(left_s.size(), right_s.size()); i++) {
            if (left_s[i] != right_s[i]) {
                diff[i] = '#';
            }
        }
        for(; i < std::max(left_s.size(), right_s.size()); i++) {
            diff[i] = '#';
        }

        std::stringstream ss;
        ss << "left:  " << left_s << "\n";
        ss << "right: " << right_s << "\n";
        ss << "diff:  " << diff << "\n";

        ASSERT_TRUE(left_s == right_s) << ss.str();
    }

}

TEST(BitCopy, test0) {
    auto const src = std::vector<uint8_t> {
        0b00000000,0b01101011,0b00000000,0b000'10100,0b10001101,0b10011101,0b00101'000,
    };
    auto dst = std::vector<uint8_t> {
        0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
    };
    auto test = [&](size_t dst_off, size_t src_off, size_t block_size) {
        copy_bits<size_t, uint8_t>(dst.data(), src.data(), dst_off, src_off, block_size, src.size());
    };

    /////////////////////////////////////////

    test(3, 9, 7);
    bit_compare(dst, std::vector<uint8_t> {
        0b00011010,0b11000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
    });

    test(0, 27, 3);
    bit_compare(dst, std::vector<uint8_t> {
        0b10111010,0b11000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
    });

    test(10, 27, 26);
    bit_compare(dst, std::vector<uint8_t> {
        0b10111010,0b11'101001,0b00011011,0b00111010,0b01010000,0b00000000,0b00000000,
//                      101001   00011011   00111010   0101
    });

    /////////////////////////////////////////

}

TEST(BitCopy, test1) {
    auto const src = std::vector<uint8_t> {
        0b11111111,0b11111111,0b11111111,0b11111111,0b11111111,0b11111111,0b11111111,
    };
    auto dst = std::vector<uint8_t> {
        0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
    };
    auto test = [&](size_t dst_off, size_t src_off, size_t block_size) {
        copy_bits<size_t, uint8_t>(dst.data(), src.data(), dst_off, src_off, block_size, src.size());
    };

    /////////////////////////////////////////

    test(3, 9, 25);
    bit_compare(dst, std::vector<uint8_t> {
        0b00011111,0b11111111,0b11111111,0b11110000,0b00000000,0b00000000,0b00000000,
    });

    test(34, 34, 16);
    bit_compare(dst, std::vector<uint8_t> {
        0b00011111,0b11111111,0b11111111,0b11110000,0b00111111,0b11111111,0b11000000,
    });

    /////////////////////////////////////////

}

TEST(BitCopy, test2) {
    auto const src = std::vector<uint8_t> {
        0b11111111,0b11111111,0b11111111,0b11111111,0b11111111,0b11111111,0b11111111,
    };
    auto dst = std::vector<uint8_t> {
        0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,0b00000000,
    };
    auto test = [&](size_t dst_off, size_t src_off, size_t block_size) {
        copy_bits<size_t, uint8_t>(dst.data(), src.data(), dst_off, src_off, block_size, src.size());
    };

    /////////////////////////////////////////

    test(32, 32, 16);
    bit_compare(dst, std::vector<uint8_t> {
        0b00000000,0b00000000,0b00000000,0b00000000,0b11111111,0b11111111,0b00000000,
    });

    /////////////////////////////////////////

}

TEST(BitCopy, test3) {
    std::cout << "size_t size: " << sizeof(size_t) * CHAR_BIT << "\n";
}
