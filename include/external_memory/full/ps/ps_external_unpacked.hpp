/*******************************************************************************
 * include/util/ps_external.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

// MATRIX WITHOUT WORD PACKING
template <typename InputType, typename stats_type>
void wx_ps_fe_builder<InputType, stats_type, false, 0>::build(
    const InputType& text,
    wavelet_structure_external& result,
    stats_type& stats) {
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

  stats.phase("l" + std::to_string(0));
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
    stats.phase("l" + std::to_string(i));
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

  stats.phase("l" + std::to_string(levels - 1));
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
template <typename InputType, typename stats_type>
void wx_ps_fe_builder<InputType, stats_type, true, 0>::build(
    const InputType& text,
    wavelet_structure_external& result,
    stats_type& stats) {
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

  result_writer_type result_writer(bvs);

  stats.phase("l0");

  PSE_VERBOSE std::cout << "Level 1 of " << levels << " (initial scan)... "
                        << std::endl;
  // Initial Scan:
  {
    leftCur->reserve(size / (64 / levels) + 1);
    rightCur->reserve(size / (64 / levels) + 1);

    input_reader_type initialReader(text);
    writer_type leftWriter(*leftCur, levels - 1);
    writer_type rightWriter(*rightCur, levels - 1);

    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned i = 0; i < 64; ++i) {
        auto symbol = *initialReader;
        auto bit = (symbol >> (levels - 1)) & 0x1;
        if (bit)
          rightWriter << symbol;
        else
          leftWriter << symbol;
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
          rightWriter << symbol;
        else
          leftWriter << symbol;
        word <<= 1;
        word |= bit;
        ++initialReader;
        hist[levels][symbol]++;
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }
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

    std::swap(leftCur, leftPrev);
    std::swap(rightCur, rightPrev);

    leftCur->clear();
    rightCur->clear();
    leftCur->reserve(size);
    rightCur->reserve(size);

    leftReader = new reader_type(*leftPrev);
    rightReader = new reader_type(*rightPrev);
    writer_type leftWriter(*leftCur);
    writer_type rightWriter(*rightCur);

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
          rightWriter << symbol;
        else
          leftWriter << symbol;
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
          rightWriter << symbol;
        else
          leftWriter << symbol;
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
  }

  stats.phase("l" + std::to_string(levels - 1));
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

//*****************************************************************************/