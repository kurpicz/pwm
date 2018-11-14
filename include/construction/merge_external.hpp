/*******************************************************************************
 * include/util/merge.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <bitset>
#include "arrays/stxxl_helper.hpp"

template <typename writer_type>
class external_merger {
  using value_type = typename writer_type::value_type;
  using vector_type = typename writer_type::vector_type;
  using reader_type = typename vector_type::bufreader_type;
  using iter_type = typename vector_type ::const_iterator;

  writer_type &writer_;
  const vector_type &vector_;
  const iter_type vector_begin_;

  uint64_t write_word = 0;
  uint8_t write_used = 0;
  uint8_t write_not_used = 64;

  std::vector<uint64_t> left_masks_vec;
  uint64_t * left_mask;

public:
  external_merger(writer_type& writer,
                  const vector_type& vector)
      : writer_(writer),
        vector_(vector),
        vector_begin_(vector_.begin()) {
    static_assert(sizeof(value_type) == 8);

    for(uint8_t width = 0; width < 64; width++) {
      left_masks_vec.push_back(((0x1ULL << width) - 1) << (64 - width));
    }
    left_masks_vec.push_back(~uint64_t(0));
    left_mask = left_masks_vec.data();
//    for(uint8_t width = 0; width < 64 + 1; width++) {
//      std::cout << std::bitset<64>(left_mask[width]).to_string() << std::endl;
//    }
//    std::cout << std::endl;
//    std::abort();
  }


  inline void writeInfix(uint64_t word,
                         const uint8_t start,
                         const uint8_t length) {
    word = (word << start) & left_mask[length];
    if(length > write_not_used) {
      writer_ << (write_word | (word >> write_used));
      write_word = word << write_not_used;
      write_used = length - write_not_used;
    } else if (length < write_not_used) {
      write_word |= (word >> write_used);
      write_used += length;
    } else {
      writer_ << (write_word | (word >> write_used));
      write_word = 0;
      write_used = 0;
    }
    write_not_used = 64 - write_used;
  }

  inline void write(const uint64_t bitcount, const uint64_t offset) {
    if(bitcount == 0) return;

    const uint64_t initial_word_id = offset / 64;
    const uint64_t final_word_id = (offset + bitcount - 1) / 64;

    const uint8_t initial_read_offset = offset % 64;
    const uint8_t initial_read_length =
        std::min(bitcount, uint64_t(64) - initial_read_offset);

    writeInfix(vector_[initial_word_id],
               initial_read_offset,
               initial_read_length);

    if(final_word_id > initial_word_id) {
      reader_type reader(vector_begin_ + initial_word_id + 1, // skip first
                         vector_begin_ + final_word_id);      // skip last
      if(write_used > 0) {
        uint64_t read_word;
        while(!reader.empty()) {
          reader >> read_word;
          writer_ << uint64_t(write_word | (read_word >> write_used));
          write_word = (read_word << write_not_used);
        }
      } else {
        for(;!reader.empty(); ++reader) writer_ << *reader;
      }

      constexpr uint8_t final_read_offset = 0;
      const uint8_t final_read_length =
          (bitcount - initial_read_length + 63) % 64 + 1;
//          bitcount - initial_read_length - 64 * (final_word_id - initial_word_id - 1);

      writeInfix(vector_[final_word_id],
                 final_read_offset,
                 final_read_length);
    }
  }

  inline void finish() {
    if(write_used > 0)
      writer_ << write_word;
    writer_.finish();
  }


};

/******************************************************************************/
