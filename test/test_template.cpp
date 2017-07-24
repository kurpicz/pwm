#include <gtest/gtest.h>
#include "test/util.hpp"

#include "benchmark/algorithm.hpp"
#include "util/common.hpp"
#include "util/debug.hpp"

// template <typename AlphabetType>
// auto ConstructWX(std::vector<AlphabetType>& text, const bool already_reduced) {
//   uint64_t max_char = 0;
//   if (!already_reduced) {
//     std::unordered_map<AlphabetType, uint64_t> word_list;
//     for (const AlphabetType& character : text) {
//     auto result = word_list.find(character);
//     if (result == word_list.end()) {
//       word_list.emplace(character, max_char++);
//     }
//     }
//     --max_char;
//     for (uint64_t i = 0; i < text.size(); ++i) {
//     text[i] = static_cast<AlphabetType>(word_list.find(text[i])->second);
//     }
//   } else {
//     for (const AlphabetType& character : text) {
//     if (character > max_char) {
//       max_char = character;
//     }
//     }
//   }
//   //std::cout << "Greatest character: " << static_cast<uint64_t>(max_char) << std::endl;

//   AlphabetType levels = 0;
//   while (max_char) {
//     max_char >>= 1;
//     ++levels;
//   }
//   //std::cout << "Levels in WM: " << static_cast<uint64_t>(levels) << std::endl;

//   return WX_NAME<AlphabetType, uint32_t>(text, text.size(), levels);
// }

TEST(wavelet_construction, smoketest) {
  auto& algo_list = algorithm_list::get_algorithm_list();
  for (const auto& a : algo_list) {
    test::roundtrip_batch([&](const std::string& s){
      auto vec = std::vector<uint8_t>(s.begin(), s.end());
      // TODO: Make computation (and reduction of text work)
      auto bvz = a->compute_bitvector(&vec, 2 , 2);
      if (a->is_tree()) {
        auto decoded_s = decode_wt(bvz.first, vec.size());
        ASSERT_EQ(s, decoded_s);
      } else {
        auto decoded_s = decode_wm(bvz.first, bvz.second, vec.size());
        ASSERT_EQ(s, decoded_s);
      }
    });
  }
}
