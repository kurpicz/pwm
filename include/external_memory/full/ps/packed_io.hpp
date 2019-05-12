/*******************************************************************************
 * include/util/ps_external.hpp
 *
 * Copyright (C) 2019 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/


#pragma once

#include "external_memory/stxxl_helper.hpp"

template <typename ValueType, bool padding>
struct word_packed_reader;

template <typename ValueType, bool padding>
struct word_packed_writer;

template <typename ValueType>
struct word_packed_reader<ValueType, true> {
  const unsigned bits;
  const uint64_t size;
  const uint64_t mask;

  const unsigned values_per_word;
  const unsigned values_in_last_word;

  stxxlvector<uint64_t>::bufreader_type reader;

  unsigned current_word_counter;
  unsigned current_shift;
  uint64_t current_word;

  word_packed_reader(const stxxlvector<uint64_t>& vec,
                     uint64_t size,
                     unsigned bits)
      : bits(bits),
        size(size),
        mask((1ULL << bits) - 1),
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
      reader >> current_word;
    }

    ValueType result = ((current_word >> current_shift) & mask);
    current_word_counter++;
    current_shift -= bits;
    return result;
  }

  inline bool empty() {
    return (reader.empty() && current_word_counter >= values_in_last_word);
  }
};

template <typename ValueType>
struct word_packed_writer<ValueType, true> {
  const unsigned bits;
  const unsigned values_per_word;
  const unsigned final_shift;
  const uint64_t mask;

  stxxlvector<uint64_t>::bufwriter_type writer;

  uint64_t total_counter;

  unsigned current_word_counter;
  uint64_t current_word;

  bool has_finished = false;

  word_packed_writer(stxxlvector<uint64_t>& vec, unsigned bits)
      : bits(bits),
        values_per_word(64 / bits),
        final_shift(64 % bits),
        mask((1ULL << bits) - 1),
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

  ~word_packed_writer() {
    finish();
  }
};

template <typename ValueType>
struct word_packed_reader<ValueType, false> {
  using value_type = ValueType;
  const int8_t bits_per_value;
  const uint64_t value_count;
  const uint64_t mask;
  const int8_t empty_bits_in_tail;

  stxxlvector<uint64_t>::bufreader_type reader;

  uint64_t current_word;
  int8_t current_shift;

  word_packed_reader(const stxxlvector<uint64_t>& vec,
                     uint64_t size,
                     uint8_t bits)
      : bits_per_value(bits),
        value_count(size),
        mask((1ULL << bits) - 1),
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
};

template <typename ValueType>
struct word_packed_writer<ValueType, false> {
  using value_type = ValueType;
  const uint8_t bits_per_value;
  const uint8_t values_per_word;
  const uint64_t mask;

  stxxlvector<uint64_t>::bufwriter_type writer;

  uint64_t word_count;
  uint64_t current_word;
  uint8_t current_shift_left;

  bool has_finished = false;

  word_packed_writer(stxxlvector<uint64_t>& vec, uint8_t bits)
      : bits_per_value(bits),
        values_per_word(64 / bits),
        mask((1ULL << bits) - 1),
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

  ~word_packed_writer() {
    finish();
  }
};


//*****************************************************************************/
