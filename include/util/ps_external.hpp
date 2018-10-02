/*******************************************************************************
 * include/util/ps_external.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

template <typename AlphabetType, typename ContextType, typename InputType>
external_bit_vectors ps_out_external(const InputType& text, uint64_t const size, const uint64_t levels,
  ContextType& ctx) {
  
  auto sorted_text_vec = std::vector<AlphabetType>(size);
  AlphabetType* const sorted_text = sorted_text_vec.data();
  
  external_bit_vectors result(levels, size);
  
  using stxxl_vector_type = typename std::remove_reference<decltype(external_bit_vectors().raw_data())>::type;
  using stxxl_writer_type = typename stxxl_vector_type::bufwriter_type;
  
  uint64_t cur_max_char = (1 << levels);

  auto& zeros = ctx.zeros();
  auto& borders = ctx.borders();
  auto& level_offsets = result.level_offsets();
  
  //stxxl vector of uint64_t
  stxxl_vector_type& bv = result.raw_data();  
  
  // While initializing the histogram, we also compute the first level
  uint64_t cur_pos = 0;
  {
    stxxl_writer_type writer(bv);
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (uint64_t i = 0; i < 64; ++i) {
        ++ctx.hist(levels, text[cur_pos + i]);
        word <<= 1;
        word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
      }
      writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (uint64_t i = 0; i < size - cur_pos; ++i) {
        ++ctx.hist(levels, text[cur_pos + i]);
        word <<= 1;
        word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
      }
      word <<= (64 - (size & 63ULL));
      writer << word;
    }
  }

  // The number of 0s at the last level is the number of "even" characters
  if (ContextType::compute_zeros) {
    for (uint64_t i = 0; i < cur_max_char; i += 2) {
      zeros[levels - 1] += ctx.hist(levels, i);
    }
  }

  // Now we compute the WM bottom-up, i.e., the last level first
  for (uint64_t level = levels - 1; level > 0; --level) {
    // Update the maximum value of a feasible a bit prefix and update the
    // histogram of the bit prefixes
    cur_max_char >>= 1;
    for (uint64_t i = 0; i < cur_max_char; ++i) {
      ctx.hist(level, i)
        = ctx.hist(level + 1, i << 1)
        + ctx.hist(level + 1, (i << 1) + 1);
    }

    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    borders[0] = 0;
    for (uint64_t i = 1; i < cur_max_char; ++i) {
      auto const prev_rho = ctx.rho(level, i - 1);

      borders[ctx.rho(level, i)] =
        borders[prev_rho] + ctx.hist(level, prev_rho);

      if (ContextType::compute_rho)  {
        ctx.set_rho(level - 1, i - 1, prev_rho >> 1);
      }
    }

    // The number of 0s is the position of the first 1 in the previous level
    if (ContextType::compute_zeros) {
      zeros[level - 1] = borders[1];
    }

    // Now we sort the text utilizing counting sort and the starting positions
    // that we have computed before
    for (uint64_t i = 0; i < size; ++i) {
      const auto cur_char = text[i];
      sorted_text[borders[cur_char >> (levels - level)]++] = cur_char;
    }

    {
      stxxl_writer_type writer(bv.begin() + level_offsets[level]);
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

  if (levels > 1) { // TODO check condition
    ctx.hist(0, 0) = ctx.hist(1, 0) + ctx.hist(1, 1);
  }


  //~ if (ContextType::compute_zeros) {
    result.setZeros(zeros);
  //~ }

  //~ std::cout << "RESULT IS: " << bv.size() << ", SHOULD BE: " << levels * ((size + 63) / 64) << ", SIZE: " << size << ", LEVELS: " << levels << std::endl << std::endl;
  //~ uint64_t words = (size + 63) / 64;
  //~ unsigned tesla = 0;
  //~ for(auto word : bv) {
    //~ if(tesla % words == 0) std::cout << std::endl << "Level " << tesla / words << ":  ";
    //~ std::cout << std::bitset<64>(word);
    //~ tesla++;
  //~ }
  //~ std::cout << std::endl << std::endl << std::endl;
  
  return result;
}


//~ template <typename ContextType, typename InputType>
//~ void ps_out_external2(const InputType& text, uint64_t const size, const uint64_t levels,
  //~ ContextType& ctx) {
  
  //~ using stxxl_vector_type = typename std::remove_reference<decltype(ctx.bv().raw_data())>::type;
  //~ using stxxl_writer_type = typename stxxl_vector_type::bufwriter_type;
  
  //~ uint64_t cur_max_char = (1 << levels);

  //~ auto& zeros = ctx.zeros();
  //~ auto& borders = ctx.borders();
  //~ auto& level_offsets = ctx.bv().level_offsets();
  
  //~ std::vector<bool> sorted_text_level_bits(size);
    
  //~ //stxxl vector of uint64_t
  //~ stxxl_vector_type& bv = ctx.bv().raw_data();  
  
  //~ // While initializing the histogram, we also compute the first level
  //~ uint64_t cur_pos = 0;
  //~ {
    //~ stxxl_writer_type writer(bv);
    //~ for (; cur_pos + 64 <= size; cur_pos += 64) {
      //~ uint64_t word = 0ULL;
      //~ for (uint64_t i = 0; i < 64; ++i) {
        //~ ++ctx.hist(levels, text[cur_pos + i]);
        //~ word <<= 1;
        //~ word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
      //~ }
      //~ writer << word;
    //~ }
    //~ if (size & 63ULL) {
      //~ uint64_t word = 0ULL;
      //~ for (uint64_t i = 0; i < size - cur_pos; ++i) {
        //~ ++ctx.hist(levels, text[cur_pos + i]);
        //~ word <<= 1;
        //~ word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
      //~ }
      //~ word <<= (64 - (size & 63ULL));
      //~ writer << word;
    //~ }
  //~ }

  //~ // The number of 0s at the last level is the number of "even" characters
  //~ if (ContextType::compute_zeros) {
    //~ for (uint64_t i = 0; i < cur_max_char; i += 2) {
      //~ zeros[levels - 1] += ctx.hist(levels, i);
    //~ }
  //~ }

  //~ // Now we compute the WM bottom-up, i.e., the last level first
  //~ for (uint64_t level = levels - 1; level > 0; --level) {
    //~ // Update the maximum value of a feasible a bit prefix and update the
    //~ // histogram of the bit prefixes
    //~ cur_max_char >>= 1;
    //~ for (uint64_t i = 0; i < cur_max_char; ++i) {
      //~ ctx.hist(level, i)
        //~ = ctx.hist(level + 1, i << 1)
        //~ + ctx.hist(level + 1, (i << 1) + 1);
    //~ }

    //~ // Compute the starting positions of characters with respect to their
    //~ // bit prefixes and the bit-reversal permutation
    //~ borders[0] = 0;
    //~ for (uint64_t i = 1; i < cur_max_char; ++i) {
      //~ auto const prev_rho = ctx.rho(level, i - 1);

      //~ borders[ctx.rho(level, i)] =
        //~ borders[prev_rho] + ctx.hist(level, prev_rho);

      //~ if (ContextType::compute_rho)  {
        //~ ctx.set_rho(level - 1, i - 1, prev_rho >> 1);
      //~ }
    //~ }

    //~ // The number of 0s is the position of the first 1 in the previous level
    //~ if (ContextType::compute_zeros) {
      //~ zeros[level - 1] = borders[1];
    //~ }

    //~ // Now we sort the text utilizing counting sort and the starting positions
    //~ // that we have computed before
    //~ for (uint64_t i = 0; i < size; ++i) {
      //~ const auto cur_char = text[i];
      //~ sorted_text_level_bits[borders[cur_char >> (levels - level)]++] 
        //~ = (cur_char >> ((levels - 1) - level)) & 1ULL;
    //~ }

    //~ {
      //~ stxxl_writer_type writer(bv.begin() + level_offsets[level]);
      //~ // Since we have sorted the text, we can simply scan it from left to right
      //~ // and for the character at position $i$ we set the $i$-th bit in the
      //~ // bit vector accordingly
      //~ for (cur_pos = 0; cur_pos + 63 < size; cur_pos += 64) {
        //~ uint64_t word = 0ULL;
        //~ for (uint64_t i = 0; i < 64; ++i) {
          //~ word <<= 1;
          //~ word |= sorted_text_level_bits[cur_pos + i];
        //~ }
        //~ writer << word;
      //~ }
      //~ if (size & 63ULL) {
        //~ uint64_t word = 0ULL;
        //~ for (uint64_t i = 0; i < size - cur_pos; ++i) {
          //~ word <<= 1;
          //~ word |= sorted_text_level_bits[cur_pos + i];
        //~ }
        //~ word <<= (64 - (size & 63ULL));
        //~ writer << word;
      //~ }
    //~ }
  //~ }

  //~ if (levels > 1) { // TODO check condition
    //~ ctx.hist(0, 0) = ctx.hist(1, 0) + ctx.hist(1, 1);
  //~ }
//~ }

//~ template <typename AlphabetType, typename ContextType, typename InputType>
//~ void ps_out_external3(const InputType& text, uint64_t const size, const uint64_t levels,
  //~ ContextType& ctx, AlphabetType* const sorted_text) {
  
  //~ using stxxl_vector_type = typename std::remove_reference<decltype(ctx.bv().raw_data())>::type;
  //~ using stxxl_writer_type = typename stxxl_vector_type::bufwriter_type;
  //~ using borders_type = typename std::remove_reference<decltype(ctx.borders())>::type;
  //~ uint64_t borders_size = ctx.borders().size();
  
  //~ uint64_t cur_max_char = (1 << levels);

  //~ auto& zeros = ctx.zeros();
  //~ auto& level_offsets = ctx.bv().level_offsets();
  //~ std::vector<borders_type> borders(levels);
  
  //~ //stxxl vector of uint64_t
  //~ stxxl_vector_type& bv = ctx.bv().raw_data();  
  
  //~ // While initializing the histogram, we also compute the first level
  //~ uint64_t cur_pos = 0;
  //~ {
    //~ stxxl_writer_type writer(bv);
    //~ for (; cur_pos + 64 <= size; cur_pos += 64) {
      //~ uint64_t word = 0ULL;
      //~ for (uint64_t i = 0; i < 64; ++i) {
        //~ ++ctx.hist(levels, text[cur_pos + i]);
        //~ word <<= 1;
        //~ word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
      //~ }
      //~ writer << word;
    //~ }
    //~ if (size & 63ULL) {
      //~ uint64_t word = 0ULL;
      //~ for (uint64_t i = 0; i < size - cur_pos; ++i) {
        //~ ++ctx.hist(levels, text[cur_pos + i]);
        //~ word <<= 1;
        //~ word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
      //~ }
      //~ word <<= (64 - (size & 63ULL));
      //~ writer << word;
    //~ }
  //~ }

  //~ // The number of 0s at the last level is the number of "even" characters
  //~ if (ContextType::compute_zeros) {
    //~ for (uint64_t i = 0; i < cur_max_char; i += 2) {
      //~ zeros[levels - 1] += ctx.hist(levels, i);
    //~ }
  //~ }

  //~ // Now we compute the WM bottom-up, i.e., the last level first
  //~ for (uint64_t level = levels - 1; level > 0; --level) {
    //~ // Update the maximum value of a feasible a bit prefix and update the
    //~ // histogram of the bit prefixes
    //~ cur_max_char >>= 1;
    //~ for (uint64_t i = 0; i < cur_max_char; ++i) {
      //~ ctx.hist(level, i)
        //~ = ctx.hist(level + 1, i << 1)
        //~ + ctx.hist(level + 1, (i << 1) + 1);
    //~ }

    //~ borders[level].resize(borders_size);
    //~ // Compute the starting positions of characters with respect to their
    //~ // bit prefixes and the bit-reversal permutation
    //~ borders[level][0] = 0;
    //~ for (uint64_t i = 1; i < cur_max_char; ++i) {
      //~ auto const prev_rho = ctx.rho(level, i - 1);

      //~ borders[level][ctx.rho(level, i)] =
        //~ borders[level][prev_rho] + ctx.hist(level, prev_rho);

      //~ if (ContextType::compute_rho)  {
        //~ ctx.set_rho(level - 1, i - 1, prev_rho >> 1);
      //~ }
    //~ }

    //~ // The number of 0s is the position of the first 1 in the previous level
    //~ if (ContextType::compute_zeros) {
      //~ zeros[level - 1] = borders[level][1];
    //~ }
  //~ }

  //~ // Now we sort the text utilizing counting sort and the starting positions
  //~ // that we have computed before
  //~ for (uint64_t i = 0; i < size; ++i) {
    //~ const auto cur_char = text[i];
    //~ for (uint64_t level = levels - 1; level > 0; --level) {
      //~ sorted_text[borders[level][cur_char >> (levels - level)]++] |= 
        //~ ((cur_char >> ((levels - 1) - level)) & 1ULL) << ((levels - 1) - level);
    //~ }
  //~ }

  //~ for (uint64_t level = levels - 1; level > 0; --level) {
    //~ {
      //~ stxxl_writer_type writer(bv.begin() + level_offsets[level]);
      //~ // Since we have sorted the text, we can simply scan it from left to right
      //~ // and for the character at position $i$ we set the $i$-th bit in the
      //~ // bit vector accordingly
      //~ for (cur_pos = 0; cur_pos + 63 < size; cur_pos += 64) {
        //~ uint64_t word = 0ULL;
        //~ for (uint64_t i = 0; i < 64; ++i) {
          //~ word <<= 1;
          //~ word |= ((sorted_text[cur_pos + i] >> ((levels - 1) - level)) & 1ULL);
        //~ }
        //~ writer << word;
      //~ }
      //~ if (size & 63ULL) {
        //~ uint64_t word = 0ULL;
        //~ for (uint64_t i = 0; i < size - cur_pos; ++i) {
          //~ word <<= 1;
          //~ word |= ((sorted_text[cur_pos + i] >> ((levels - 1) - level)) & 1ULL);
        //~ }
        //~ word <<= (64 - (size & 63ULL));
        //~ writer << word;
      //~ }
    //~ }
  //~ }

  //~ if (levels > 1) { // TODO check condition
    //~ ctx.hist(0, 0) = ctx.hist(1, 0) + ctx.hist(1, 1);
  //~ }
//~ }


template <typename ReaderType, typename WriterType>
inline void split(
    ReaderType &reader, 
    WriterType &leftWriter, 
    WriterType &rightWriter,
    unsigned shift, 
    uint64_t length) {
  for(uint64_t i = 0; i < length; i++) {
      auto symbol = *reader;
      if((symbol >> shift) & 0x1)
        rightWriter << symbol;
      else
        leftWriter << symbol;
      ++reader;
  }
}

// DO NOT USE
template <typename AlphabetType, typename ContextType, typename InputType>
external_bit_vectors ps_fully_external_matrix2(const InputType& text, uint64_t const size, const uint64_t levels,
  ContextType& /*ctx*/) {
  
  std::cout << "PS external" << std::endl;
  external_bit_vectors result(levels, size, 0);
  std::vector<uint64_t>& zeros = result.zeros();
  zeros.resize(levels);

  using in_vector_type = InputType;
  using reader_type = typename in_vector_type::bufreader_type;
  using writer_type = typename in_vector_type::bufwriter_type;
  
  using out_vector_type = typename std::remove_reference<decltype(result.raw_data())>::type;
  using result_writer_type = typename out_vector_type::bufwriter_type;
  out_vector_type& bv = result.raw_data();

  std::vector<std::vector<uint64_t>> hist(levels + 1);
  for(unsigned i = 0; i <= levels; i++)
    hist[i].resize(pow(2, i));

  stxxl_files::reset_usage();
  in_vector_type v1 = stxxl_files::getVector<in_vector_type>(1);
  in_vector_type v2 = stxxl_files::getVector<in_vector_type>(2);
  in_vector_type v3 = stxxl_files::getVector<in_vector_type>(3);
  in_vector_type v4 = stxxl_files::getVector<in_vector_type>(4);

  v1.clear();
  v2.clear();
  v3.clear();
  v4.clear();
  
  in_vector_type * leftPrev = &v1;
  in_vector_type * rightPrev = &v2;
  
  in_vector_type * leftCur = &v3;
  in_vector_type * rightCur = &v4;

  leftCur->reserve(size);
  rightCur->reserve(size);
  
  reader_type * leftReader;
  reader_type * rightReader;
  writer_type * leftWriter;
  writer_type * rightWriter;
    
  result_writer_type result_writer(bv);
  
  std::cout << "Level 1 of " << levels << " (initial scan)... " << std::endl;
  // Initial Scan:
  {
    leftReader = new reader_type(text);
    leftWriter = new writer_type(*leftCur);
    rightWriter = new writer_type(*rightCur);
    
    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned i = 0; i < 64; ++i) {
        auto symbol = **leftReader;
        if(symbol >> (levels - 1))
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        ++(*leftReader);
        hist[levels][symbol]++;
        word <<= 1;
        word |= ((symbol >> (levels - 1)) & 1ULL);
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned i = 0; i < size - cur_pos; ++i) {
        auto symbol = **leftReader;
        if(symbol >> (levels - 1))
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        ++(*leftReader);
        hist[levels][symbol]++;
        word <<= 1;
        word |= ((symbol >> (levels - 1)) & 1ULL);
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    delete leftReader;
    delete leftWriter;
    delete rightWriter;
  }
  
  // calculate histograms
  for(unsigned i = levels; i > 0; i--) {
    for(unsigned j = 0; j < hist[i].size(); j+=2) {
      hist[i - 1][j / 2] = hist[i][j] + hist[i][j + 1];
    }
  }
  
  // scans (top down WT construction in left-right-buffers)
  for(unsigned i = 1; i < levels; i++) {
    std::cout << "Level " << i + 1 << " of " << levels << "... " << std::endl;
    if(i > 1) {
      leftCur->reserve(size);
      rightCur->reserve(size);

      leftReader = new reader_type(*leftPrev);
      rightReader = new reader_type(*rightPrev);
      leftWriter = new writer_type(*leftCur);
      rightWriter = new writer_type(*rightCur);
    
      for(unsigned j = 0; j < hist[i - 1].size(); j+=2) {
        split(*leftReader, *leftWriter, *rightWriter, levels - i, hist[i - 1][j]);
        split(*rightReader, *leftWriter, *rightWriter, levels - i, hist[i - 1][j + 1]);
      }
    
      delete leftReader;
      delete rightReader;
      delete leftWriter;
      delete rightWriter;
    }
    zeros[i - 1] = leftCur->size();
    std::cout << "ZEROES " << zeros[i-1] << "  " << size << std::endl;

    std::swap(leftCur, leftPrev);
    std::swap(rightCur, rightPrev);
    
    leftCur->clear();
    rightCur->clear();
    
    leftReader = new reader_type(*leftPrev);
    rightReader = new reader_type(*rightPrev);

    reader_type * leftRightReader = leftReader;
    
    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < 64; ++k) {
        if(leftRightReader->empty()) {
          leftRightReader = rightReader;
        }
        auto symbol = **leftRightReader;
        word <<= 1;
        word |= ((symbol >> (levels - i - 1)) & 1ULL);
        
        ++(*leftRightReader);
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < size - cur_pos; ++k) {
        if(leftRightReader->empty()) {
          leftRightReader = rightReader;
        }
        auto symbol = **leftRightReader;
        word <<= 1;
        word |= ((symbol >> (levels - i - 1)) & 1ULL);
        
        ++(*leftRightReader);
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }
    
    delete leftReader;
    delete rightReader;
    
  }
  
//  delete leftCur;
//  delete rightCur;
//  delete leftPrev;
//  delete rightPrev;
  
  result_writer.finish();
  
  std::cout << "Done." << std::endl << std::endl;
  
  //~ std::cout << "RESULT IS: " << bv.size() << ", SHOULD BE: " << levels * ((size + 63) / 64) << ", SIZE: " << size << ", LEVELS: " << levels << std::endl << std::endl;
  //~ uint64_t words = (size + 63) / 64;
  //~ unsigned tesla = 0;
  //~ for(auto word : bv) {
    //~ if(tesla % words == 0) std::cout << std::endl << "Level " << tesla / words << ":  ";
    //~ std::cout << std::bitset<64>(word);
    //~ tesla++;
  //~ }
  //~ std::cout << std::endl << std::endl << std::endl;
  
  return result;
}


//template <typename ReaderType, typename WriterType>
//inline void split(
//    ReaderType &reader,
//    WriterType &leftWriter,
//    WriterType &rightWriter,
//    unsigned shift,
//    uint64_t length) {
//  for(uint64_t i = 0; i < length; i++) {
//    auto symbol = *reader;
//    if((symbol >> shift) & 0x1)
//      rightWriter << symbol;
//    else
//      leftWriter << symbol;
//    ++reader;
//  }
//}


template <typename AlphabetType, typename ContextType, typename InputType>
external_bit_vectors ps_fully_external_matrix(const InputType& text, uint64_t const size, const uint64_t levels,
                                            ContextType& /*ctx*/) {

  std::cout << "PS external" << std::endl;
  external_bit_vectors result(levels, size, 0);
  std::vector<uint64_t>& zeros = result.zeros();
  zeros.resize(levels);

  using in_vector_type = InputType;
  using reader_type = typename in_vector_type::bufreader_type;
  using writer_type = typename in_vector_type::bufwriter_type;

  using out_vector_type = typename std::remove_reference<decltype(result.raw_data())>::type;
  using result_writer_type = typename out_vector_type::bufwriter_type;
  out_vector_type& bv = result.raw_data();

  stxxl_files::reset_usage();
  in_vector_type v1 = stxxl_files::getVector<in_vector_type>(1);
  in_vector_type v2 = stxxl_files::getVector<in_vector_type>(2);
  in_vector_type v3 = stxxl_files::getVector<in_vector_type>(3);
  in_vector_type v4 = stxxl_files::getVector<in_vector_type>(4);

  v1.clear();
  v2.clear();
  v3.clear();
  v4.clear();

  in_vector_type * leftPrev = &v1;
  in_vector_type * rightPrev = &v2;

  in_vector_type * leftCur = &v3;
  in_vector_type * rightCur = &v4;

  reader_type * leftReader;
  reader_type * rightReader;
  writer_type * leftWriter;
  writer_type * rightWriter;

  result_writer_type result_writer(bv);

  // scans (top down WT construction in left-right-buffers)
  for(unsigned i = 0; i < levels; i++) {
    std::cout << "Level " << i + 1 << " of " << levels << "... " << std::endl;

    std::swap(leftCur, leftPrev);
    std::swap(rightCur, rightPrev);

    leftCur->clear();
    rightCur->clear();
    leftCur->reserve(size);
    rightCur->reserve(size);

    if(i == 0) {
      leftReader = new reader_type(text);
    } else {
      zeros[i - 1] = leftPrev->size();
      leftReader = new reader_type(*leftPrev);
    }

    rightReader = new reader_type(*rightPrev);
    leftWriter = new writer_type(*leftCur);
    rightWriter = new writer_type(*rightCur);

    reader_type * leftRightReader = leftReader;

    auto const shift = levels - i - 1;
    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < 64; ++k) {
        if(leftRightReader->empty()) {
          leftRightReader = rightReader;
        }
        auto symbol = **leftRightReader;
        if((symbol >> shift) & 0x1)
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        word <<= 1;
        word |= ((symbol >> shift) & 1ULL);

        ++(*leftRightReader);
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < size - cur_pos; ++k) {
        if(leftRightReader->empty()) {
          leftRightReader = rightReader;
        }
        auto symbol = **leftRightReader;
        if((symbol >> shift) & 0x1)
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        word <<= 1;
        word |= ((symbol >> shift) & 1ULL);

        ++(*leftRightReader);
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    delete leftReader;
    delete rightReader;
    delete leftWriter;
    delete rightWriter;
  }

  result_writer.finish();

  std::cout << "Done." << std::endl << std::endl;

  //~ std::cout << "RESULT IS: " << bv.size() << ", SHOULD BE: " << levels * ((size + 63) / 64) << ", SIZE: " << size << ", LEVELS: " << levels << std::endl << std::endl;
  //~ uint64_t words = (size + 63) / 64;
  //~ unsigned tesla = 0;
  //~ for(auto word : bv) {
  //~ if(tesla % words == 0) std::cout << std::endl << "Level " << tesla / words << ":  ";
  //~ std::cout << std::bitset<64>(word);
  //~ tesla++;
  //~ }
  //~ std::cout << std::endl << std::endl << std::endl;

  return result;
}


template <typename AlphabetType, typename ContextType, typename InputType>
external_bit_vectors ps_fully_external_tree(const InputType& text, uint64_t const size, const uint64_t levels,
                                       ContextType& /*ctx*/) {

  std::cout << "PS external" << std::endl;
  external_bit_vectors result(levels, size, 0);

  using in_vector_type = InputType;
  using reader_type = typename in_vector_type::bufreader_type;
  using writer_type = typename in_vector_type::bufwriter_type;

  using out_vector_type = typename std::remove_reference<decltype(result.raw_data())>::type;
  using result_writer_type = typename out_vector_type::bufwriter_type;
  out_vector_type& bv = result.raw_data();

  std::vector<std::vector<uint64_t>> hist(levels + 1);
  for(unsigned i = 0; i <= levels; i++)
    hist[i].resize(pow(2, i));

  stxxl_files::reset_usage();
  in_vector_type v1 = stxxl_files::getVector<in_vector_type>(1);
  in_vector_type v2 = stxxl_files::getVector<in_vector_type>(2);
  in_vector_type v3 = stxxl_files::getVector<in_vector_type>(3);
  in_vector_type v4 = stxxl_files::getVector<in_vector_type>(4);

  v1.clear();
  v2.clear();
  v3.clear();
  v4.clear();

  in_vector_type * leftPrev = &v1;
  in_vector_type * rightPrev = &v2;

  in_vector_type * leftCur = &v3;
  in_vector_type * rightCur = &v4;

  reader_type * leftReader;
  reader_type * rightReader;
  writer_type * leftWriter;
  writer_type * rightWriter;

  result_writer_type result_writer(bv);

  std::cout << "Level 1 of " << levels << " (initial scan)... " << std::endl;
  // Initial Scan:
  {
    leftCur->reserve(size);
    rightCur->reserve(size);

    leftReader = new reader_type(text);
    leftWriter = new writer_type(*leftCur);
    rightWriter = new writer_type(*rightCur);

    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned i = 0; i < 64; ++i) {
        auto symbol = **leftReader;
        if(symbol >> (levels - 1))
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        ++(*leftReader);
        hist[levels][symbol]++;
        word <<= 1;
        word |= ((symbol >> (levels - 1)) & 1ULL);
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned i = 0; i < size - cur_pos; ++i) {
        auto symbol = **leftReader;
        if(symbol >> (levels - 1))
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        ++(*leftReader);
        hist[levels][symbol]++;
        word <<= 1;
        word |= ((symbol >> (levels - 1)) & 1ULL);
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    delete leftReader;
    delete leftWriter;
    delete rightWriter;
  }

  // calculate histograms
  for(unsigned i = levels; i > 0; i--) {
    for(unsigned j = 0; j < hist[i].size(); j+=2) {
      hist[i - 1][j / 2] = hist[i][j] + hist[i][j + 1];
    }
  }

  // scans (top down WT construction in left-right-buffers)
  for(unsigned i = 1; i < levels; i++) {
    std::cout << "Level " << i + 1 << " of " << levels << "... " << std::endl;

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

    unsigned histId = 0;
    uint64_t histRemains = hist[i][0];
    reader_type * leftRightReader = leftReader;

    auto const shift = levels - i - 1;
    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < 64; ++k) {
        if(histRemains == 0) {
          do{
            histId++;
            histRemains = hist[i][histId];
          } while(histRemains == 0);
          if(histId % 2 == 0) leftRightReader = leftReader;
          else leftRightReader = rightReader;
        }
        auto symbol = **leftRightReader;
        if((symbol >> shift) & 0x1)
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        word <<= 1;
        word |= ((symbol >> shift) & 1ULL);

        ++(*leftRightReader);
        --histRemains;
      }
      result_writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < size - cur_pos; ++k) {
        if(histRemains == 0) {
          do{
            histId++;
            histRemains = hist[i][histId];
          } while(histRemains == 0);
          if(histId % 2 == 0) leftRightReader = leftReader;
          else leftRightReader = rightReader;
        }
        auto symbol = **leftRightReader;
        if((symbol >> shift) & 0x1)
          *rightWriter << symbol;
        else
          *leftWriter << symbol;
        word <<= 1;
        word |= ((symbol >> shift) & 1ULL);

        ++(*leftRightReader);
        --histRemains;
      }
      word <<= (64 - (size & 63ULL));
      result_writer << word;
    }

    delete leftReader;
    delete rightReader;
    delete leftWriter;
    delete rightWriter;
  }

  result_writer.finish();

  std::cout << "Done." << std::endl << std::endl;

  //~ std::cout << "RESULT IS: " << bv.size() << ", SHOULD BE: " << levels * ((size + 63) / 64) << ", SIZE: " << size << ", LEVELS: " << levels << std::endl << std::endl;
  //~ uint64_t words = (size + 63) / 64;
  //~ unsigned tesla = 0;
  //~ for(auto word : bv) {
  //~ if(tesla % words == 0) std::cout << std::endl << "Level " << tesla / words << ":  ";
  //~ std::cout << std::bitset<64>(word);
  //~ tesla++;
  //~ }
  //~ std::cout << std::endl << std::endl << std::endl;

  return result;
}



/******************************************************************************/
