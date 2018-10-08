/*******************************************************************************
 * include/util/ps_external.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <thread>
#include <mutex>
#include <condition_variable>

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




struct {
  uint64_t get(uint8_t width) {
    uint64_t mask = 0ULL;
    for(unsigned i = 0; i < width; i++) {
      mask <<= 1;
      mask |= 1ULL;
    }
    return mask;
  }
} mask_factory;


template <typename ValueType>
struct word_packed_reader_na {
  using value_type = ValueType;
  const int8_t bits_per_value;
  const uint64_t value_count;
  const uint64_t mask;
  const int8_t empty_bits_in_tail;

  stxxlvector<uint64_t>::bufreader_type reader;

  uint64_t current_word;
  int8_t current_shift;

  word_packed_reader_na(const stxxlvector<uint64_t>& vec, uint64_t size, uint8_t bits) :
      bits_per_value(bits),
      value_count(size),
      mask(mask_factory.get(bits)),
      empty_bits_in_tail(64 - ((size * bits + 63) % 64 + 1)),
      reader(vec),
      current_word(0),
      current_shift(-bits) {

    assert((size * bits_per_value) / 64 <= vec.size());
    assert((size * bits_per_value) / 64 + 1 >= size);
  }

  ValueType next() {
    if(current_shift < 0) {
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

  uint64_t size() const {
    return value_count;
  }

  bool empty() const {
    return reader.empty() && current_shift < empty_bits_in_tail;
  }

  ~word_packed_reader_na() {

  }
};


template <typename ValueType>
struct word_packed_writer_na {
  const uint8_t bits_per_value;
  const uint8_t values_per_word;
  const uint64_t mask;

  stxxlvector<uint64_t>::bufwriter_type writer;

  uint64_t word_count;
  uint64_t current_word;
  uint8_t current_shift_left;

  bool has_finished = false;

  word_packed_writer_na(stxxlvector<uint64_t>& vec, uint8_t bits) :
      bits_per_value(bits),
      values_per_word(64 / bits),
      mask(mask_factory.get(bits)),
      writer(vec),
      word_count(1),
      current_word(0),
      current_shift_left(64) {
  }

  void next(ValueType value) {
    if(current_shift_left < bits_per_value) {
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

  uint64_t size() const {
    return (word_count * 64 - current_shift_left) / bits_per_value;
  }

  uint64_t finish() {
    if(has_finished) return size();
    has_finished = true;
    if (current_shift_left != 64) {
      current_word <<= current_shift_left;
      writer << current_word;
      writer.finish();
    }
    return size();
  }

  ~word_packed_writer_na() {
    finish();
  }

};


template <typename reader_type>
struct buffered_wp_reader {
  
  using value_type = typename reader_type::value_type;
  reader_type& reader;

  const uint64_t value_count;
  const uint64_t values_per_buffer;
  const uint64_t tail_count;
  uint64_t swaps_left;
  uint64_t it;

  value_type * front_buffer;
  value_type * back_buffer;

  buffered_wp_reader(reader_type &reader, uint64_t bytes) :
      reader(reader),
      value_count(reader.size()),
      values_per_buffer((bytes / sizeof(value_type)) / 2),
      tail_count(value_count == 0 ? 0 :
        (value_count + values_per_buffer - 1) % values_per_buffer + 1),
      swaps_left((value_count + values_per_buffer - 1) / values_per_buffer),
      it(values_per_buffer) {
    front_buffer = new value_type [values_per_buffer];
    back_buffer = new value_type [values_per_buffer];
    if(swaps_left > 0) load_back_buffer();
  }

  void load_back_buffer() {
    const uint64_t count = (swaps_left == 1) ? tail_count : values_per_buffer;
    for(uint64_t i = 0; i < count; ++i) {
      back_buffer[i] = reader.next();
    }
  }

  value_type next() {
    if(it < values_per_buffer) {
      return front_buffer[it++];
    } else {
      --swaps_left;
      std::swap(front_buffer, back_buffer);
      if(swaps_left > 0) load_back_buffer();
      it = 1;
      return front_buffer[0];
    }
  }

  bool empty() const {
    return swaps_left == 0 && it >= tail_count;
  }

  uint64_t size() const {
    return value_count;
  }

  ~buffered_wp_reader() {
    delete [] front_buffer;
    delete [] back_buffer;
  }

};


template <typename ValueType>
struct word_packed_reader {
  const unsigned bits;
  const uint64_t size;

  const unsigned values_per_word;
  const unsigned values_in_last_word;

  stxxlvector<uint64_t>::bufreader_type reader;

  unsigned current_word_counter;
  unsigned current_shift;
  uint64_t current_word;
  uint64_t mask;

  word_packed_reader(const stxxlvector<uint64_t>& vec, uint64_t size, unsigned bits) :
      bits(bits),
      size(size),
      values_per_word(64 / bits),
      values_in_last_word((size - 1) % (64 / bits) + 1),
      reader(vec) {

    assert((vec.size() - 1) * values_per_word < size);
    assert(vec.size() * values_per_word >= size);

    current_word_counter = 65;

    mask = 0ULL;
    for(unsigned i = 0; i < bits; i++) {
      mask <<= 1;
      mask |= 1ULL;
    }
  }

  ValueType next() {
    if(current_word_counter >= values_per_word) {
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

  bool empty() {
    return
      (reader.empty() &&
       current_word_counter >= values_in_last_word);
  }

  ~word_packed_reader() {

  }
};


template <typename ValueType>
struct word_packed_writer {
  const unsigned bits;
  const unsigned values_per_word;
  const unsigned final_shift;

  stxxlvector<uint64_t>::bufwriter_type writer;

  uint64_t total_counter;

  unsigned current_word_counter;
  uint64_t current_word;
  uint64_t mask;

  bool has_finished = false;

  word_packed_writer(stxxlvector<uint64_t>& vec, unsigned bits) :
      bits(bits),
      values_per_word(64 / bits),
      final_shift(64 % bits),
      writer(vec) {

    total_counter = 0;
    current_word_counter = 0;
    current_word = 0ULL;

    mask = 0ULL;
    for(unsigned i = 0; i < bits; i++) {
      mask <<= 1;
      mask |= 1ULL;
    }
  }

  void next(ValueType value) {
    current_word <<= bits;
    current_word |= (mask & value);
    current_word_counter++;

    if(current_word_counter == values_per_word) {
      total_counter += values_per_word;
      writer << (current_word << final_shift);
      current_word_counter = 0;
      current_word = 0ULL;
    }
  }

  uint64_t size() const {
    return total_counter;
  }

  uint64_t finish() {
    if(has_finished) return size();
    auto final_size = total_counter + current_word_counter;
    while(current_word_counter > 0)
      next(ValueType(0));
    total_counter = final_size;
    writer.finish();
    return size();
  }

  ~word_packed_writer() {
    finish();
  }

};


template <typename AlphabetType, typename ContextType, typename InputType>
external_bit_vectors ps_fully_external_matrix(const InputType& text, uint64_t const size, const uint64_t levels,
                                            ContextType& /*ctx*/) {

  std::cout << "PS external" << std::endl;
  external_bit_vectors result(levels, size, 0);
  std::vector<uint64_t>& zeros = result.zeros();
  zeros.resize(levels);

  uint64_t bytes_per_reader = 1024 * 1024 * 32;

  using input_reader_type = typename InputType::bufreader_type;
  using out_vector_type = typename std::remove_reference<decltype(result.raw_data())>::type;
  using result_writer_type = typename out_vector_type::bufwriter_type;

  using word_packed_vector_type = stxxlvector<uint64_t>;
  using value_type = typename InputType::value_type;
  using reader_type = word_packed_reader_na<value_type>;
  using writer_type = word_packed_writer_na<value_type>;
  using bufreader_type = buffered_wp_reader<reader_type>;

  out_vector_type& bv = result.raw_data();

  stxxl_files::reset_usage();
  word_packed_vector_type v1 = stxxl_files::getVector<word_packed_vector_type>(1);
  word_packed_vector_type v2 = stxxl_files::getVector<word_packed_vector_type>(2);
  word_packed_vector_type v3 = stxxl_files::getVector<word_packed_vector_type>(3);
  word_packed_vector_type v4 = stxxl_files::getVector<word_packed_vector_type>(4);

  v1.clear();
  v2.clear();
  v3.clear();
  v4.clear();

  word_packed_vector_type * leftPrev = &v1;
  word_packed_vector_type * rightPrev = &v2;

  word_packed_vector_type * leftCur = &v3;
  word_packed_vector_type * rightCur = &v4;

  bufreader_type * leftReader;
  bufreader_type * rightReader;
  writer_type * leftWriter;
  writer_type * rightWriter;

  result_writer_type result_writer(bv);

  uint64_t left_size = 0;
  uint64_t right_size = 0;

  std::cout << "Level 1 of " << levels << " (initial scan)... " << std::endl;
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
        if(bit)
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
        if(bit)
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
  for(unsigned i = 1; i < levels - 1; i++) {
    std::cout << "Level " << i + 1 << " of " << levels << "... " << std::endl;

    const unsigned bits = levels - i;
    const unsigned values_per_word = 64 / bits;

    std::swap(leftCur, leftPrev);
    std::swap(rightCur, rightPrev);

    leftCur->clear();
    rightCur->clear();
    leftCur->reserve(size / values_per_word + 1);
    rightCur->reserve(size / values_per_word + 1);

    zeros[i - 1] = left_size;
    reader_type lRead(*leftPrev, left_size, bits);
    reader_type rRead(*rightPrev, right_size, bits);
    leftReader = new bufreader_type(lRead, bytes_per_reader);
    rightReader = new bufreader_type(rRead, bytes_per_reader);
    leftWriter = new writer_type(*leftCur, bits - 1);
    rightWriter = new writer_type(*rightCur, bits - 1);

    //buffered_wp_reader<reader_type> testreader(*leftReader, 16);

    bufreader_type * leftRightReader = leftReader;

    auto const shift = levels - i - 1;
    uint64_t cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (unsigned k = 0; k < 64; ++k) {
        if(leftRightReader->empty()) {
          leftRightReader = rightReader;
        }
        auto symbol = leftRightReader->next();
        auto bit = (symbol >> shift) & 0x1;
        if(bit)
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
        if(leftRightReader->empty()) {
          leftRightReader = rightReader;
        }
        auto symbol = leftRightReader->next();
        auto bit = (symbol >> shift) & 0x1;
        if(bit)
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

  std::cout << "Level " << levels << " of " << levels << " (final scan)... " << std::endl;

  zeros[levels - 2] = left_size;
  reader_type lRead(*leftCur, left_size, 1);
  reader_type rRead(*rightCur, right_size, 1);
  leftReader = new bufreader_type(lRead, bytes_per_reader);
  rightReader = new bufreader_type(rRead, bytes_per_reader);

  bufreader_type * leftRightReader = leftReader;

  uint64_t cur_pos = 0;
  for (; cur_pos + 64 <= size; cur_pos += 64) {
    uint64_t word = 0ULL;
    for (unsigned k = 0; k < 64; ++k) {
      if(leftRightReader->empty()) {
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
      if(leftRightReader->empty()) {
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

  std::cout << "Done." << std::endl << std::endl;

  return result;
}


template <typename AlphabetType, typename ContextType, typename InputType>
external_bit_vectors ps_fully_external_tree(const InputType& text, uint64_t const size, const uint64_t levels,
                                       ContextType& /*ctx*/) {

  std::cout << "PS external" << std::endl;
  external_bit_vectors result(levels, size, 0);

  using input_reader_type = typename InputType::bufreader_type;
  using out_vector_type = typename std::remove_reference<decltype(result.raw_data())>::type;
  using result_writer_type = typename out_vector_type::bufwriter_type;

  using word_packed_vector_type = stxxlvector<uint64_t>;
  using value_type = typename InputType::value_type;
  using reader_type = word_packed_reader_na<value_type>;
  using writer_type = word_packed_writer_na<value_type>;

  out_vector_type& bv = result.raw_data();

  stxxl_files::reset_usage();
  word_packed_vector_type v1 = stxxl_files::getVector<word_packed_vector_type>(1);
  word_packed_vector_type v2 = stxxl_files::getVector<word_packed_vector_type>(2);
  word_packed_vector_type v3 = stxxl_files::getVector<word_packed_vector_type>(3);
  word_packed_vector_type v4 = stxxl_files::getVector<word_packed_vector_type>(4);

  v1.clear();
  v2.clear();
  v3.clear();
  v4.clear();

  word_packed_vector_type * leftPrev = &v1;
  word_packed_vector_type * rightPrev = &v2;

  word_packed_vector_type * leftCur = &v3;
  word_packed_vector_type * rightCur = &v4;

  reader_type * leftReader;
  reader_type * rightReader;
  writer_type * leftWriter;
  writer_type * rightWriter;

  result_writer_type result_writer(bv);

  uint64_t left_size = 0;
  uint64_t right_size = 0;

  std::vector<std::vector<uint64_t>> hist(levels + 1);
  for(unsigned i = 0; i <= levels; i++)
    hist[i].resize(pow(2, i));

  std::cout << "Level 1 of " << levels << " (initial scan)... " << std::endl;
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
        if(bit)
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
        if(bit)
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
  for(unsigned i = levels; i > 0; i--) {
    for(unsigned j = 0; j < hist[i].size(); j+=2) {
      hist[i - 1][j / 2] = hist[i][j] + hist[i][j + 1];
    }
  }

  // scans (top down WT construction in left-right-buffers)
  for(unsigned i = 1; i < levels - 1; i++) {
    std::cout << "Level " << i + 1 << " of " << levels << "... " << std::endl;

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

    reader_type * leftRightReader = leftReader;

    unsigned histId = 0;
    uint64_t histRemains = hist[i][0];
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
        auto symbol = leftRightReader->next();
        auto bit = (symbol >> shift) & 0x1;
        if(bit)
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
        if(histRemains == 0) {
          do{
            histId++;
            histRemains = hist[i][histId];
          } while(histRemains == 0);
          if(histId % 2 == 0) leftRightReader = leftReader;
          else leftRightReader = rightReader;
        }
        auto symbol = leftRightReader->next();
        auto bit = (symbol >> shift) & 0x1;
        if(bit)
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

  std::cout << "Level " << levels << " of " << levels << " (final scan)... " << std::endl;

  leftReader = new reader_type(*leftCur, left_size, 1);
  rightReader = new reader_type(*rightCur, right_size, 1);

  reader_type * leftRightReader = leftReader;

  unsigned histId = 0;
  uint64_t histRemains = hist[levels - 1][0];
  uint64_t cur_pos = 0;
  for (; cur_pos + 64 <= size; cur_pos += 64) {
    uint64_t word = 0ULL;
    for (unsigned k = 0; k < 64; ++k) {
      if(histRemains == 0) {
        do{
          histId++;
          histRemains = hist[levels - 1][histId];
        } while(histRemains == 0);
        if(histId % 2 == 0) leftRightReader = leftReader;
        else leftRightReader = rightReader;
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
      if(histRemains == 0) {
        do{
          histId++;
          histRemains = hist[levels - 1][histId];
        } while(histRemains == 0);
        if(histId % 2 == 0) leftRightReader = leftReader;
        else leftRightReader = rightReader;
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

  std::cout << "Done." << std::endl << std::endl;

  return result;
}



/******************************************************************************/
