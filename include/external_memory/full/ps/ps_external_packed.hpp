/*******************************************************************************
 * include/util/ps_external.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

// MATRIX WITH WORD PACKING
template <typename InputType, typename stats_type, int word_packing_mode>
void wx_ps_fe_builder<InputType, stats_type, false, word_packing_mode>::build(
    const InputType& text,
    wavelet_structure_external& result,
    stats_type& stats) {
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
  using reader_type = word_packed_reader<value_type, (word_packing_mode == 1)>;
  using writer_type = word_packed_writer<value_type, (word_packing_mode == 1)>;

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

  result_writer_type result_writer(bvs);

  uint64_t left_size = 0;
  uint64_t right_size = 0;

  stats.phase("l" + std::to_string(0));
  PSE_VERBOSE std::cout << "Level 1 of " << levels << " (initial scan)... "
                        << std::endl;
  {
    leftCur->reserve(size / (64 / levels) + 1);
    rightCur->reserve(size / (64 / levels) + 1);

    input_reader_type initialReader(text);
    writer_type leftWriter(*leftCur, levels - 1);
    writer_type rightWriter(*rightCur, levels - 1);

    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < 64; ++k) {
        auto symbol = *initialReader;
        auto bit = (symbol >> (levels - 1)) & 0x1;
        if (bit)
          rightWriter.next(symbol);
        else
          leftWriter.next(symbol);
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
          rightWriter.next(symbol);
        else
          leftWriter.next(symbol);
        word <<= 1;
        word |= bit;
        ++initialReader;
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    left_size = leftWriter.finish();
    right_size = rightWriter.finish();
  }

  // scans (top down WT construction in left-right-buffers)
  for (unsigned i = 1; i < levels - 1; i++) {
    stats.phase("l" + std::to_string(i));
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
    writer_type leftWriter(*leftCur, bits - 1);
    writer_type rightWriter(*rightCur, bits - 1);

    // buffered_wp_reader<reader_type> testreader(*leftReader, 16);

    reader_type* leftRightReader = leftReader;

    auto const shift = levels - i - 1;
    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < 64; ++k) {
        if (PWM_UNLIKELY(leftRightReader->empty())) {
          leftRightReader = rightReader;
        }
        auto symbol = leftRightReader->next();
        auto bit = (symbol >> shift) & 0x1;
        if (bit)
          rightWriter.next(symbol);
        else
          leftWriter.next(symbol);
        word <<= 1;
        word |= bit;
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < size - cur_pos; ++k) {
        if (PWM_UNLIKELY(leftRightReader->empty())) {
          leftRightReader = rightReader;
        }
        auto symbol = leftRightReader->next();
        auto bit = (symbol >> shift) & 0x1;
        if (bit)
          rightWriter.next(symbol);
        else
          leftWriter.next(symbol);
        word <<= 1;
        word |= bit;
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    left_size = leftWriter.finish();
    right_size = rightWriter.finish();
    delete leftReader;
    delete rightReader;
  }

  stats.phase("l" + std::to_string(levels - 1));
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
      if (PWM_UNLIKELY(leftRightReader->empty())) {
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
      if (PWM_UNLIKELY(leftRightReader->empty())) {
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
template <typename InputType, typename stats_type, int word_packing_mode>
void wx_ps_fe_builder<InputType, stats_type, true, word_packing_mode>::build(
    const InputType& text,
    wavelet_structure_external& result,
    stats_type& stats) {
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
  using reader_type = word_packed_reader<value_type, (word_packing_mode == 1)>;
  using writer_type = word_packed_writer<value_type, (word_packing_mode == 1)>;

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

  stats.phase("l" + std::to_string(0));
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
    stats.phase("l" + std::to_string(i));
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

  stats.phase("l" + std::to_string(levels - 1));
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

//*****************************************************************************/