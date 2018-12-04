/*******************************************************************************
 * include/util/ps_external.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "construction/ctx_single_level_external.hpp"
#include "construction/wavelet_structure_external.hpp"
#include "util/debug_assert.hpp"

#define PSE_VERBOSE if (false)

template <typename AlphabetType, bool is_tree, typename InputType>
void ps_out_external(
    const InputType& text,
    wavelet_structure_external& result) {
  auto levels = result.levels();
  auto size = result.text_size();

  auto& bvs = wavelet_structure_external_writer::bvs(result);
  auto& zeros = wavelet_structure_external_writer::zeros(result);
  auto& hist = wavelet_structure_external_writer::histograms(result);
  auto const& level_offsets = result.level_offsets();

  ctx_single_level_external<is_tree> ctx(size, levels);
  auto& borders = ctx.borders();

  auto sorted_text_vec = std::vector<AlphabetType>(size);
  AlphabetType* sorted_text = sorted_text_vec.data();

  uint64_t cur_max_char = (1 << levels);

  using result_writer_type =
      typename std::remove_reference<decltype(bvs)>::type::bufwriter_type;

  // While initializing the histogram, we also compute the first level
  uint64_t cur_pos = 0;
  {
    result_writer_type writer(bvs);
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (uint64_t i = 0; i < 64; ++i) {
        ++hist[levels][text[cur_pos + i]];
        word <<= 1;
        word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
      }
      writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (uint64_t i = 0; i < size - cur_pos; ++i) {
        ++hist[levels][text[cur_pos + i]];
        word <<= 1;
        word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
      }
      word <<= (64 - (size & 63ULL));
      writer << word;
    }
  }

  // The number of 0s at the last level is the number of "even" characters
  if constexpr (!is_tree) {
    for (uint64_t i = 0; i < cur_max_char; i += 2) {
      zeros[levels - 1] += hist[levels][i];
    }
  }

  // Now we compute the WM bottom-up, i.e., the last level first
  for (uint64_t level = levels - 1; level > 0; --level) {
    // Update the maximum value of a feasible a bit prefix and update the
    // histogram of the bit prefixes
    cur_max_char >>= 1;
    for (uint64_t i = 0; i < cur_max_char; ++i) {
      hist[level][i] =
          hist[level + 1][i << 1] + hist[level + 1][(i << 1) + 1];
    }

    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    borders[0] = 0;
    for (uint64_t i = 1; i < cur_max_char; ++i) {
      auto const prev_rho = ctx.rho(level, i - 1);

      borders[ctx.rho(level, i)] =
          borders[prev_rho] + hist[level][prev_rho];

      if constexpr (!is_tree) {
        ctx.set_rho(level - 1, i - 1, prev_rho >> 1);
      }
    }

    // The number of 0s is the position of the first 1 in the previous level
    if constexpr (!is_tree) {
      zeros[level - 1] = borders[1];
    }

    // Now we sort the text utilizing counting sort and the starting positions
    // that we have computed before
    for (uint64_t i = 0; i < size; ++i) {
      const auto cur_char = text[i];
      sorted_text[borders[cur_char >> (levels - level)]++] = cur_char;
    }

    {
      result_writer_type writer(bvs.begin() + level_offsets[level]);
      // Since we have sorted the text, we can simply scan it from left to right
      // and for the character at position $i$ we set the $i$-th bit in the
      // bit vector accordingly
      for (cur_pos = 0; cur_pos + 63 < size; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < 64; ++i) {
          word <<= 1;
          word |= ((sorted_text[cur_pos + i] >> ((levels - 1) - level)) & 1ULL);
        }
        writer << word;
      }
      if (size & 63ULL) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < size - cur_pos; ++i) {
          word <<= 1;
          word |= ((sorted_text[cur_pos + i] >> ((levels - 1) - level)) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        writer << word;
      }
    }
  }
}

struct {
  uint64_t get(uint8_t width) {
    uint64_t mask = 0ULL;
    for (unsigned i = 0; i < width; i++) {
      mask <<= 1;
      mask |= 1ULL;
    }
    return mask;
  }
} mask_factory;

template <typename ValueType>
struct word_packed_reader_with_padding {
  const unsigned bits;
  const uint64_t size;
  const uint64_t mask;

  const unsigned values_per_word;
  const unsigned values_in_last_word;

  stxxlvector<uint64_t>::bufreader_type reader;

  unsigned current_word_counter;
  unsigned current_shift;
  uint64_t current_word;

  word_packed_reader_with_padding(const stxxlvector<uint64_t>& vec,
                                  uint64_t size,
                                  unsigned bits)
      : bits(bits),
        size(size),
        mask(mask_factory.get(bits)),
        values_per_word(64 / bits),
        values_in_last_word((size - 1) % (64 / bits) + 1),
        reader(vec) {
    DCHECK(vec.size() * values_per_word < size + values_per_word);
    DCHECK(vec.size() * values_per_word >= size);

    current_word_counter = 65;
  }

  inline ValueType next() {
    if (current_word_counter >= values_per_word) {
      current_word_counter = 0;
      current_shift = 64 - bits;
      current_word = *reader;
      ++reader;
    }

    ValueType result = ((current_word >> current_shift) & mask);
    current_word_counter++;
    current_shift -= bits;
    return result;
  }

  inline bool empty() {
    return (reader.empty() && current_word_counter >= values_in_last_word);
  }

  ~word_packed_reader_with_padding() {}
};

template <typename ValueType>
struct word_packed_writer_with_padding {
  const unsigned bits;
  const unsigned values_per_word;
  const unsigned final_shift;
  const uint64_t mask;

  stxxlvector<uint64_t>::bufwriter_type writer;

  uint64_t total_counter;

  unsigned current_word_counter;
  uint64_t current_word;

  bool has_finished = false;

  word_packed_writer_with_padding(stxxlvector<uint64_t>& vec, unsigned bits)
      : bits(bits),
        values_per_word(64 / bits),
        final_shift(64 % bits),
        mask(mask_factory.get(bits)),
        writer(vec) {

    total_counter = 0;
    current_word_counter = 0;
    current_word = 0ULL;
  }

  inline void next(ValueType value) {
    current_word <<= bits;
    current_word |= (mask & value);
    current_word_counter++;

    if (current_word_counter == values_per_word) {
      total_counter += values_per_word;
      writer << (current_word << final_shift);
      current_word_counter = 0;
      current_word = 0ULL;
    }
  }

  inline uint64_t size() const {
    return total_counter;
  }

  inline uint64_t finish() {
    if (has_finished)
      return size();
    auto final_size = total_counter + current_word_counter;
    while (current_word_counter > 0)
      next(ValueType(0));
    total_counter = final_size;
    writer.finish();
    return size();
  }

  ~word_packed_writer_with_padding() {
    finish();
  }
};

template <typename ValueType>
struct word_packed_reader_without_padding {
  using value_type = ValueType;
  const int8_t bits_per_value;
  const uint64_t value_count;
  const uint64_t mask;
  const int8_t empty_bits_in_tail;

  stxxlvector<uint64_t>::bufreader_type reader;

  uint64_t current_word;
  int8_t current_shift;

  word_packed_reader_without_padding(const stxxlvector<uint64_t>& vec,
                                     uint64_t size,
                                     uint8_t bits)
      : bits_per_value(bits),
        value_count(size),
        mask(mask_factory.get(bits)),
        empty_bits_in_tail(64 - ((size * bits + 63) % 64 + 1)),
        reader(vec),
        current_word(0),
        current_shift(-bits) {

    DCHECK((size * bits_per_value) / 64 <= vec.size());
    DCHECK((size * bits_per_value) / 64 + 1 >= vec.size());
  }

  inline ValueType next() {
    if (current_shift < 0) {
      ValueType result = (current_word << -current_shift) & mask;
      current_word = *reader;
      current_shift += 64;
      result |= current_word >> current_shift;
      current_shift -= bits_per_value;
      ++reader;
      return result;
    }

    ValueType result = ((current_word >> current_shift) & mask);
    current_shift -= bits_per_value;
    return result;
  }

  inline uint64_t size() const {
    return value_count;
  }

  inline bool empty() const {
    return reader.empty() && current_shift < empty_bits_in_tail;
  }

  ~word_packed_reader_without_padding() {}
};

template <typename ValueType>
struct word_packed_writer_without_padding {
  using value_type = ValueType;
  const uint8_t bits_per_value;
  const uint8_t values_per_word;
  const uint64_t mask;

  stxxlvector<uint64_t>::bufwriter_type writer;

  uint64_t word_count;
  uint64_t current_word;
  uint8_t current_shift_left;

  bool has_finished = false;

  word_packed_writer_without_padding(stxxlvector<uint64_t>& vec, uint8_t bits)
      : bits_per_value(bits),
        values_per_word(64 / bits),
        mask(mask_factory.get(bits)),
        writer(vec),
        word_count(1),
        current_word(0),
        current_shift_left(64) {}

  inline void next(ValueType value) {
    if (current_shift_left < bits_per_value) {
      current_word <<= current_shift_left;
      current_shift_left = bits_per_value - current_shift_left;
      current_word |= ((value & mask) >> current_shift_left);
      current_shift_left = 64 - current_shift_left;
      writer << current_word;
      current_word = value;
      ++word_count;
      return;
    }

    current_word <<= bits_per_value;
    current_word |= (mask & value);
    current_shift_left -= bits_per_value;
  }

  inline uint64_t size() const {
    return (word_count * 64 - current_shift_left) / bits_per_value;
  }

  inline uint64_t finish() {
    if (has_finished)
      return size();
    has_finished = true;
    if (current_shift_left != 64) {
      current_word <<= current_shift_left;
      writer << current_word;
      writer.finish();
    }
    return size();
  }

  ~word_packed_writer_without_padding() {
    finish();
  }
};

template <typename value_type, int word_packing_mode>
struct reader_writer_types {};

template <typename value_type>
struct reader_writer_types<value_type, 1> {
  using reader_type = word_packed_reader_with_padding<value_type>;
  using writer_type = word_packed_writer_with_padding<value_type>;
};

template <typename value_type>
struct reader_writer_types<value_type, 2> {
  using reader_type = word_packed_reader_without_padding<value_type>;
  using writer_type = word_packed_writer_without_padding<value_type>;
};

template <typename InputType, bool is_tree, int word_packing_mode = 2>
struct wx_ps_fe_builder {
  static void build(const InputType& text, wavelet_structure_external& result);
};

template <typename InputType, int word_packing_mode>
struct wx_ps_fe_builder<InputType, false, word_packing_mode> {
  static void build(const InputType& text, wavelet_structure_external& result);
};

template <typename InputType, int word_packing_mode>
struct wx_ps_fe_builder<InputType, true, word_packing_mode> {
  static void build(const InputType& text, wavelet_structure_external& result);
};

template <typename InputType>
struct wx_ps_fe_builder<InputType, false, 0> {
  static void build(const InputType& text, wavelet_structure_external& result);
};

template <typename InputType>
struct wx_ps_fe_builder<InputType, true, 0> {
  static void build(const InputType& text, wavelet_structure_external& result);
};

// MATRIX WITH WORD PACKING
template <typename InputType, int word_packing_mode>
void wx_ps_fe_builder<InputType, false, word_packing_mode>::build(
    const InputType& text, wavelet_structure_external& result) {
  PSE_VERBOSE std::cout << "PS external (MATRIX, WORDPACKING "
                        << word_packing_mode << ")" << std::endl;
  auto levels = result.levels();
  auto size = result.text_size();

  auto& zeros = wavelet_structure_external_writer::zeros(result);
  auto& bvs = wavelet_structure_external_writer::bvs(result);

  using input_reader_type = typename InputType::bufreader_type;
  using result_writer_type =
      typename std::remove_reference<decltype(bvs)>::type::bufwriter_type;

  using value_type = typename InputType::value_type;
  using word_packed_vector_type = stxxlvector<uint64_t>;
  using reader_type = word_packed_reader_with_padding<value_type>;
  using writer_type = word_packed_writer_with_padding<value_type>;

  word_packed_vector_type v1 =
      stxxl_files::getVectorTemporary<word_packed_vector_type>(1);
  word_packed_vector_type v2 =
      stxxl_files::getVectorTemporary<word_packed_vector_type>(2);
  word_packed_vector_type v3 =
      stxxl_files::getVectorTemporary<word_packed_vector_type>(3);
  word_packed_vector_type v4 =
      stxxl_files::getVectorTemporary<word_packed_vector_type>(4);

  v1.clear();
  v2.clear();
  v3.clear();
  v4.clear();

  word_packed_vector_type* leftPrev = &v1;
  word_packed_vector_type* rightPrev = &v2;

  word_packed_vector_type* leftCur = &v3;
  word_packed_vector_type* rightCur = &v4;

  reader_type* leftReader;
  reader_type* rightReader;
  writer_type* leftWriter;
  writer_type* rightWriter;

  result_writer_type result_writer(bvs);

  uint64_t left_size = 0;
  uint64_t right_size = 0;

  PSE_VERBOSE std::cout << "Level 1 of " << levels << " (initial scan)... "
                        << std::endl;
  {
    leftCur->reserve(size / (64 / levels) + 1);
    rightCur->reserve(size / (64 / levels) + 1);

    input_reader_type initialReader(text);
    leftWriter = new writer_type(*leftCur, levels - 1);
    rightWriter = new writer_type(*rightCur, levels - 1);

    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < 64; ++k) {
        auto symbol = *initialReader;
        auto bit = (symbol >> (levels - 1)) & 0x1;
        if (bit)
          rightWriter->next(symbol);
        else
          leftWriter->next(symbol);
        word <<= 1;
        word |= bit;
        ++initialReader;
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < size - cur_pos; ++k) {
        auto symbol = *initialReader;
        auto bit = (symbol >> (levels - 1)) & 0x1;
        if (bit)
          rightWriter->next(symbol);
        else
          leftWriter->next(symbol);
        word <<= 1;
        word |= bit;
        ++initialReader;
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    left_size = leftWriter->finish();
    right_size = rightWriter->finish();
    delete leftWriter;
    delete rightWriter;
  }

  // scans (top down WT construction in left-right-buffers)
  for (unsigned i = 1; i < levels - 1; i++) {
    PSE_VERBOSE std::cout << "Level " << i + 1 << " of " << levels << "... "
                          << std::endl;

    const unsigned bits = levels - i;
    const unsigned values_per_word = 64 / bits;

    std::swap(leftCur, leftPrev);
    std::swap(rightCur, rightPrev);

    leftCur->clear();
    rightCur->clear();
    leftCur->reserve(size / values_per_word + 1);
    rightCur->reserve(size / values_per_word + 1);

    zeros[i - 1] = left_size;
    leftReader = new reader_type(*leftPrev, left_size, bits);
    rightReader = new reader_type(*rightPrev, right_size, bits);
    leftWriter = new writer_type(*leftCur, bits - 1);
    rightWriter = new writer_type(*rightCur, bits - 1);

    // buffered_wp_reader<reader_type> testreader(*leftReader, 16);

    reader_type* leftRightReader = leftReader;

    auto const shift = levels - i - 1;
    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < 64; ++k) {
        if (leftRightReader->empty()) {
          leftRightReader = rightReader;
        }
        auto symbol = leftRightReader->next();
        auto bit = (symbol >> shift) & 0x1;
        if (bit)
          rightWriter->next(symbol);
        else
          leftWriter->next(symbol);
        word <<= 1;
        word |= bit;
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < size - cur_pos; ++k) {
        if (leftRightReader->empty()) {
          leftRightReader = rightReader;
        }
        auto symbol = leftRightReader->next();
        auto bit = (symbol >> shift) & 0x1;
        if (bit)
          rightWriter->next(symbol);
        else
          leftWriter->next(symbol);
        word <<= 1;
        word |= bit;
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    left_size = leftWriter->finish();
    right_size = rightWriter->finish();
    delete leftReader;
    delete rightReader;
    delete leftWriter;
    delete rightWriter;
  }

  PSE_VERBOSE std::cout << "Level " << levels << " of " << levels
                        << " (final scan)... " << std::endl;

  zeros[levels - 2] = left_size;
  leftReader = new reader_type(*leftCur, left_size, 1);
  rightReader = new reader_type(*rightCur, right_size, 1);

  reader_type* leftRightReader = leftReader;

  uint64_t cur_pos = 0;
  for (; cur_pos + 64 <= size; cur_pos += 64) {
    uint64_t word = 0ULL;
    for (unsigned k = 0; k < 64; ++k) {
      if (leftRightReader->empty()) {
        leftRightReader = rightReader;
      }
      word <<= 1;
      word |= (leftRightReader->next() & 1ULL);
    }
    result_writer << word;
  }
  if (size & 63ULL) {
    uint64_t word = 0ULL;
    for (unsigned k = 0; k < size - cur_pos; ++k) {
      if (leftRightReader->empty()) {
        leftRightReader = rightReader;
      }
      word <<= 1;
      word |= (leftRightReader->next() & 1ULL);
    }
    word <<= (64 - (size & 63ULL));
    result_writer << word;
  }

  delete leftReader;
  delete rightReader;

  result_writer.finish();

  PSE_VERBOSE std::cout << "Done." << std::endl << std::endl;
}

// TREE WITH WORD PACKING
template <typename InputType, int word_packing_mode>
void wx_ps_fe_builder<InputType, true, word_packing_mode>::build(
    const InputType& text, wavelet_structure_external& result) {
  PSE_VERBOSE std::cout << "PS external (TREE, WORDPACKING "
                        << word_packing_mode << ")" << std::endl;
  auto levels = result.levels();
  auto size = result.text_size();

  auto& hist = wavelet_structure_external_writer::histograms(result);
  auto& bvs = wavelet_structure_external_writer::bvs(result);

  using input_reader_type = typename InputType::bufreader_type;
  using result_writer_type =
      typename std::remove_reference<decltype(bvs)>::type::bufwriter_type;

  using value_type = typename InputType::value_type;
  using word_packed_vector_type = stxxlvector<uint64_t>;
  using reader_type = word_packed_reader_with_padding<value_type>;
  using writer_type = word_packed_writer_with_padding<value_type>;

  word_packed_vector_type v1 =
      stxxl_files::getVectorTemporary<word_packed_vector_type>(1);
  word_packed_vector_type v2 =
      stxxl_files::getVectorTemporary<word_packed_vector_type>(2);
  word_packed_vector_type v3 =
      stxxl_files::getVectorTemporary<word_packed_vector_type>(3);
  word_packed_vector_type v4 =
      stxxl_files::getVectorTemporary<word_packed_vector_type>(4);

  v1.clear();
  v2.clear();
  v3.clear();
  v4.clear();

  word_packed_vector_type* leftPrev = &v1;
  word_packed_vector_type* rightPrev = &v2;

  word_packed_vector_type* leftCur = &v3;
  word_packed_vector_type* rightCur = &v4;

  reader_type* leftReader;
  reader_type* rightReader;
  writer_type* leftWriter;
  writer_type* rightWriter;

  result_writer_type result_writer(bvs);

  uint64_t left_size = 0;
  uint64_t right_size = 0;

  PSE_VERBOSE std::cout << "Level 1 of " << levels << " (initial scan)... "
                        << std::endl;
  // Initial Scan:
  {
    leftCur->reserve(size / (64 / levels) + 1);
    rightCur->reserve(size / (64 / levels) + 1);

    input_reader_type initialReader(text);
    leftWriter = new writer_type(*leftCur, levels - 1);
    rightWriter = new writer_type(*rightCur, levels - 1);

    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned i = 0; i < 64; ++i) {
        auto symbol = *initialReader;
        auto bit = (symbol >> (levels - 1)) & 0x1;
        if (bit)
          rightWriter->next(symbol);
        else
          leftWriter->next(symbol);
        word <<= 1;
        word |= bit;
        ++initialReader;
        hist[levels][symbol]++;
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned i = 0; i < size - cur_pos; ++i) {
        auto symbol = *initialReader;
        auto bit = (symbol >> (levels - 1)) & 0x1;
        if (bit)
          rightWriter->next(symbol);
        else
          leftWriter->next(symbol);
        word <<= 1;
        word |= bit;
        ++initialReader;
        hist[levels][symbol]++;
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    left_size = leftWriter->finish();
    right_size = rightWriter->finish();
    delete leftWriter;
    delete rightWriter;
  }

  // calculate histograms
  for (unsigned i = levels; i > 0; i--) {
    for (unsigned j = 0; j < hist[i].size(); j += 2) {
      hist[i - 1][j / 2] = hist[i][j] + hist[i][j + 1];
    }
  }

  // scans (top down WT construction in left-right-buffers)
  for (unsigned i = 1; i < levels - 1; i++) {
    PSE_VERBOSE std::cout << "Level " << i + 1 << " of " << levels << "... "
                          << std::endl;

    const unsigned bits = levels - i;
    const unsigned values_per_word = 64 / bits;

    std::swap(leftCur, leftPrev);
    std::swap(rightCur, rightPrev);

    leftCur->clear();
    rightCur->clear();
    leftCur->reserve(size / values_per_word + 1);
    rightCur->reserve(size / values_per_word + 1);

    leftReader = new reader_type(*leftPrev, left_size, bits);
    rightReader = new reader_type(*rightPrev, right_size, bits);
    leftWriter = new writer_type(*leftCur, bits - 1);
    rightWriter = new writer_type(*rightCur, bits - 1);

    reader_type* leftRightReader = leftReader;

    unsigned histId = 0;
    uint64_t histRemains = hist[i][0];
    auto const shift = levels - i - 1;
    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < 64; ++k) {
        if (histRemains == 0) {
          do {
            histId++;
            histRemains = hist[i][histId];
          } while (histRemains == 0);
          if (histId % 2 == 0)
            leftRightReader = leftReader;
          else
            leftRightReader = rightReader;
        }
        auto symbol = leftRightReader->next();
        auto bit = (symbol >> shift) & 0x1;
        if (bit)
          rightWriter->next(symbol);
        else
          leftWriter->next(symbol);
        word <<= 1;
        word |= bit;
        --histRemains;
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < size - cur_pos; ++k) {
        if (histRemains == 0) {
          do {
            histId++;
            histRemains = hist[i][histId];
          } while (histRemains == 0);
          if (histId % 2 == 0)
            leftRightReader = leftReader;
          else
            leftRightReader = rightReader;
        }
        auto symbol = leftRightReader->next();
        auto bit = (symbol >> shift) & 0x1;
        if (bit)
          rightWriter->next(symbol);
        else
          leftWriter->next(symbol);
        word <<= 1;
        word |= bit;
        --histRemains;
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    left_size = leftWriter->finish();
    right_size = rightWriter->finish();
    delete leftReader;
    delete rightReader;
    delete leftWriter;
    delete rightWriter;
  }

  PSE_VERBOSE std::cout << "Level " << levels << " of " << levels
                        << " (final scan)... " << std::endl;

  leftReader = new reader_type(*leftCur, left_size, 1);
  rightReader = new reader_type(*rightCur, right_size, 1);

  reader_type* leftRightReader = leftReader;

  unsigned histId = 0;
  uint64_t histRemains = hist[levels - 1][0];
  uint64_t cur_pos = 0;
  for (; cur_pos + 64 <= size; cur_pos += 64) {
    uint64_t word = 0ULL;
    for (unsigned k = 0; k < 64; ++k) {
      if (histRemains == 0) {
        do {
          histId++;
          histRemains = hist[levels - 1][histId];
        } while (histRemains == 0);
        if (histId % 2 == 0)
          leftRightReader = leftReader;
        else
          leftRightReader = rightReader;
      }
      word <<= 1;
      word |= (leftRightReader->next() & 1ULL);
      --histRemains;
    }
    result_writer << word;
  }
  if (size & 63ULL) {
    uint64_t word = 0ULL;
    for (unsigned k = 0; k < size - cur_pos; ++k) {
      if (histRemains == 0) {
        do {
          histId++;
          histRemains = hist[levels - 1][histId];
        } while (histRemains == 0);
        if (histId % 2 == 0)
          leftRightReader = leftReader;
        else
          leftRightReader = rightReader;
      }
      word <<= 1;
      word |= (leftRightReader->next() & 1ULL);
      --histRemains;
    }
    word <<= (64 - (size & 63ULL));
    result_writer << word;
  }

  delete leftReader;
  delete rightReader;

  result_writer.finish();

  PSE_VERBOSE std::cout << "Done." << std::endl << std::endl;
}

// MATRIX WITHOUT WORD PACKING
template <typename InputType>
void wx_ps_fe_builder<InputType, false, 0>::build(
    const InputType& text, wavelet_structure_external& result) {
  PSE_VERBOSE std::cout << "PS external (MATRIX, NO WORDPACKING)" << std::endl;
  auto levels = result.levels();
  auto size = result.text_size();

  auto& zeros = wavelet_structure_external_writer::zeros(result);
  auto& bvs = wavelet_structure_external_writer::bvs(result);

  using input_reader_type = typename InputType::bufreader_type;
  using result_writer_type =
  typename std::remove_reference<decltype(bvs)>::type::bufwriter_type;

  using value_type = typename InputType::value_type;
  using vector_type = stxxlvector<value_type>;
  using reader_type = typename vector_type::bufreader_type;
  using writer_type = typename vector_type::bufwriter_type;

  vector_type v1 = stxxl_files::getVectorTemporary<vector_type>(1);
  vector_type v2 = stxxl_files::getVectorTemporary<vector_type>(2);
  vector_type v3 = stxxl_files::getVectorTemporary<vector_type>(3);
  vector_type v4 = stxxl_files::getVectorTemporary<vector_type>(4);

  v1.clear();
  v2.clear();
  v3.clear();
  v4.clear();

  vector_type* leftPrev = &v1;
  vector_type* rightPrev = &v2;

  vector_type* leftCur = &v3;
  vector_type* rightCur = &v4;

  reader_type* leftReader;
  reader_type* rightReader;
  writer_type* leftWriter;
  writer_type* rightWriter;

  result_writer_type result_writer(bvs);

  PSE_VERBOSE std::cout << "Level 1 of " << levels << " (initial scan)... "
                        << std::endl;
  {
    leftCur->reserve(size);
    rightCur->reserve(size);

    input_reader_type initialReader(text);
    leftWriter = new writer_type(*leftCur);
    rightWriter = new writer_type(*rightCur);

    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < 64; ++k) {
        auto symbol = *initialReader;
        auto bit = (symbol >> (levels - 1)) & 0x1;
        if (bit)
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        word <<= 1;
        word |= bit;
        ++initialReader;
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < size - cur_pos; ++k) {
        auto symbol = *initialReader;
        auto bit = (symbol >> (levels - 1)) & 0x1;
        if (bit)
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        word <<= 1;
        word |= bit;
        ++initialReader;
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    delete leftWriter;
    delete rightWriter;
  }

  // scans (top down WT construction in left-right-buffers)
  for (unsigned i = 1; i < levels - 1; i++) {
    PSE_VERBOSE std::cout << "Level " << i + 1 << " of " << levels << "... "
                          << std::endl;

    std::swap(leftCur, leftPrev);
    std::swap(rightCur, rightPrev);

    leftCur->clear();
    rightCur->clear();
    leftCur->reserve(size);
    rightCur->reserve(size);

    zeros[i - 1] = leftPrev->size();
    leftReader = new reader_type(*leftPrev);
    rightReader = new reader_type(*rightPrev);
    leftWriter = new writer_type(*leftCur);
    rightWriter = new writer_type(*rightCur);

    // buffered_wp_reader<reader_type> testreader(*leftReader, 16);

    reader_type* leftRightReader = leftReader;

    auto const shift = levels - i - 1;
    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < 64; ++k) {
        if (leftRightReader->empty()) {
          leftRightReader = rightReader;
        }
        auto symbol = **leftRightReader;
        auto bit = (symbol >> shift) & 0x1;
        if (bit)
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        word <<= 1;
        word |= bit;
        ++*leftRightReader;
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < size - cur_pos; ++k) {
        if (leftRightReader->empty()) {
          leftRightReader = rightReader;
        }
        auto symbol = **leftRightReader;
        auto bit = (symbol >> shift) & 0x1;
        if (bit)
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        word <<= 1;
        word |= bit;
        ++*leftRightReader;
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    delete leftReader;
    delete rightReader;
    delete leftWriter;
    delete rightWriter;
  }

  PSE_VERBOSE std::cout << "Level " << levels << " of " << levels
                        << " (final scan)... " << std::endl;

  zeros[levels - 2] = leftCur->size();
  leftReader = new reader_type(*leftCur);
  rightReader = new reader_type(*rightCur);

  reader_type* leftRightReader = leftReader;

  uint64_t cur_pos = 0;
  for (; cur_pos + 64 <= size; cur_pos += 64) {
    uint64_t word = 0ULL;
    for (unsigned k = 0; k < 64; ++k) {
      if (leftRightReader->empty()) {
        leftRightReader = rightReader;
      }
      word <<= 1;
      word |= **leftRightReader & 1ULL;
      ++*leftRightReader;
    }
    result_writer << word;
  }
  if (size & 63ULL) {
    uint64_t word = 0ULL;
    for (unsigned k = 0; k < size - cur_pos; ++k) {
      if (leftRightReader->empty()) {
        leftRightReader = rightReader;
      }
      word <<= 1;
      word |= **leftRightReader & 1ULL;
      ++*leftRightReader;
    }
    word <<= (64 - (size & 63ULL));
    result_writer << word;
  }

  delete leftReader;
  delete rightReader;

  result_writer.finish();

  PSE_VERBOSE std::cout << "Done." << std::endl << std::endl;
}

// TREE WITHOUT WORD PACKING
template <typename InputType>
void wx_ps_fe_builder<InputType, true, 0>::build(
    const InputType& text, wavelet_structure_external& result) {
  PSE_VERBOSE std::cout << "PS external (TREE, NO WORDPACKING)" << std::endl;
  auto levels = result.levels();
  auto size = result.text_size();

  auto& hist = wavelet_structure_external_writer::histograms(result);
  auto& bvs = wavelet_structure_external_writer::bvs(result);

  using input_reader_type = typename InputType::bufreader_type;
  using result_writer_type =
  typename std::remove_reference<decltype(bvs)>::type::bufwriter_type;

  using value_type = typename InputType::value_type;
  using vector_type = stxxlvector<value_type>;
  using reader_type = typename vector_type::bufreader_type;
  using writer_type = typename vector_type::bufwriter_type;

  vector_type v1 = stxxl_files::getVectorTemporary<vector_type>(1);
  vector_type v2 = stxxl_files::getVectorTemporary<vector_type>(2);
  vector_type v3 = stxxl_files::getVectorTemporary<vector_type>(3);
  vector_type v4 = stxxl_files::getVectorTemporary<vector_type>(4);

  v1.clear();
  v2.clear();
  v3.clear();
  v4.clear();

  vector_type* leftPrev = &v1;
  vector_type* rightPrev = &v2;

  vector_type* leftCur = &v3;
  vector_type* rightCur = &v4;

  reader_type* leftReader;
  reader_type* rightReader;
  writer_type* leftWriter;
  writer_type* rightWriter;

  result_writer_type result_writer(bvs);

  PSE_VERBOSE std::cout << "Level 1 of " << levels << " (initial scan)... "
                        << std::endl;
  // Initial Scan:
  {
    leftCur->reserve(size / (64 / levels) + 1);
    rightCur->reserve(size / (64 / levels) + 1);

    input_reader_type initialReader(text);
    leftWriter = new writer_type(*leftCur, levels - 1);
    rightWriter = new writer_type(*rightCur, levels - 1);

    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned i = 0; i < 64; ++i) {
        auto symbol = *initialReader;
        auto bit = (symbol >> (levels - 1)) & 0x1;
        if (bit)
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        word <<= 1;
        word |= bit;
        ++initialReader;
        hist[levels][symbol]++;
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned i = 0; i < size - cur_pos; ++i) {
        auto symbol = *initialReader;
        auto bit = (symbol >> (levels - 1)) & 0x1;
        if (bit)
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        word <<= 1;
        word |= bit;
        ++initialReader;
        hist[levels][symbol]++;
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    delete leftWriter;
    delete rightWriter;
  }

  // calculate histograms
  for (unsigned i = levels; i > 0; i--) {
    for (unsigned j = 0; j < hist[i].size(); j += 2) {
      hist[i - 1][j / 2] = hist[i][j] + hist[i][j + 1];
    }
  }

  // scans (top down WT construction in left-right-buffers)
  for (unsigned i = 1; i < levels - 1; i++) {
    PSE_VERBOSE std::cout << "Level " << i + 1 << " of " << levels << "... "
                          << std::endl;

    std::swap(leftCur, leftPrev);
    std::swap(rightCur, rightPrev);

    leftCur->clear();
    rightCur->clear();
    leftCur->reserve(size);
    rightCur->reserve(size);

    leftReader = new reader_type(*leftPrev);
    rightReader = new reader_type(*rightPrev);
    leftWriter = new writer_type(*leftCur);
    rightWriter = new writer_type(*rightCur);

    reader_type* leftRightReader = leftReader;

    unsigned histId = 0;
    uint64_t histRemains = hist[i][0];
    auto const shift = levels - i - 1;
    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < 64; ++k) {
        if (histRemains == 0) {
          do {
            histId++;
            histRemains = hist[i][histId];
          } while (histRemains == 0);
          if (histId % 2 == 0)
            leftRightReader = leftReader;
          else
            leftRightReader = rightReader;
        }
        auto symbol = **leftRightReader;
        auto bit = (symbol >> shift) & 0x1;
        if (bit)
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        word <<= 1;
        word |= bit;
        --histRemains;
        ++*leftRightReader;
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < size - cur_pos; ++k) {
        if (histRemains == 0) {
          do {
            histId++;
            histRemains = hist[i][histId];
          } while (histRemains == 0);
          if (histId % 2 == 0)
            leftRightReader = leftReader;
          else
            leftRightReader = rightReader;
        }
        auto symbol = **leftRightReader;
        auto bit = (symbol >> shift) & 0x1;
        if (bit)
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        word <<= 1;
        word |= bit;
        --histRemains;
        ++*leftRightReader;
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    delete leftReader;
    delete rightReader;
    delete leftWriter;
    delete rightWriter;
  }

  PSE_VERBOSE std::cout << "Level " << levels << " of " << levels
                        << " (final scan)... " << std::endl;

  leftReader = new reader_type(*leftCur);
  rightReader = new reader_type(*rightCur);

  reader_type* leftRightReader = leftReader;

  unsigned histId = 0;
  uint64_t histRemains = hist[levels - 1][0];
  uint64_t cur_pos = 0;
  for (; cur_pos + 64 <= size; cur_pos += 64) {
    uint64_t word = 0ULL;
    for (unsigned k = 0; k < 64; ++k) {
      if (histRemains == 0) {
        do {
          histId++;
          histRemains = hist[levels - 1][histId];
        } while (histRemains == 0);
        if (histId % 2 == 0)
          leftRightReader = leftReader;
        else
          leftRightReader = rightReader;
      }
      word <<= 1;
      word |= **leftRightReader & 1ULL;
      --histRemains;
      ++*leftRightReader;
    }
    result_writer << word;
  }
  if (size & 63ULL) {
    uint64_t word = 0ULL;
    for (unsigned k = 0; k < size - cur_pos; ++k) {
      if (histRemains == 0) {
        do {
          histId++;
          histRemains = hist[levels - 1][histId];
        } while (histRemains == 0);
        if (histId % 2 == 0)
          leftRightReader = leftReader;
        else
          leftRightReader = rightReader;
      }
      word <<= 1;
      word |= **leftRightReader & 1ULL;
      --histRemains;
      ++*leftRightReader;
    }
    word <<= (64 - (size & 63ULL));
    result_writer << word;
  }

  delete leftReader;
  delete rightReader;

  result_writer.finish();

  PSE_VERBOSE std::cout << "Done." << std::endl << std::endl;
}

/******************************************************************************/
