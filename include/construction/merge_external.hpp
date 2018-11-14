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

template <typename vector_type>
class external_merger {
  using value_type = typename vector_type::value_type;
  using writer_type = typename vector_type::bufwriter_type;
  using reader_type = typename vector_type::bufreader_type;
  using iter_type = typename vector_type ::const_iterator;

  writer_type writer_;
  const vector_type * vector_;
  const iter_type vector_begin_;
  const iter_type vector_end_;

  uint64_t write_word = 0;
  uint8_t write_used = 0;
  uint8_t write_not_used = 64;

  std::vector<uint64_t> left_masks_vec;
  uint64_t * left_mask;

  inline void wstr(uint64_t) {
    //std::cout << "Write: " << std::bitset<64>(word).to_string() << std::endl;
  }

  inline void rstr(uint64_t) {
    //std::cout << "Read:  " << std::bitset<64>(word).to_string() << std::endl;
  }

public:
  external_merger(const vector_type * input, vector_type& output)
      : writer_(output),
        vector_(input),
        vector_begin_(input->begin()),
        vector_end_(input->end()){
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

    //std::cout << "TEMP BVS IN WRITER (VEC): " << std::endl;
//    for(auto word : *vector_) {
      //std::cout << std::bitset<64>(word).to_string() << std::endl;
//    }

    //std::cout << "TEMP BVS IN WRITER (IT): " << std::endl;
//    for(auto it = vector_begin_; it != vector_end_; it++) {
      //std::cout << std::bitset<64>(*it).to_string() << std::endl;
//    }

    //std::cout << "ADDRESS IN WRITER: " << &output << std::endl;
  }

  inline void finishLevel() {
    if(write_used > 0) {
      wstr(write_word);
      writer_ << write_word;
      write_word = 0;
      write_used = 0;
      write_not_used = 64;
    }
  }

  inline void writeInfix(uint64_t word,
                         const uint8_t start,
                         const uint8_t length) {
    word = (word << start) & left_mask[length];
    if(length > write_not_used) {
      wstr((write_word | (word >> write_used)));
      writer_ << (write_word | (word >> write_used));
      write_word = word << write_not_used;
      write_used = length - write_not_used;
    } else if (length < write_not_used) {
      write_word |= (word >> write_used);
      write_used += length;
    } else {
      wstr((write_word | (word >> write_used)));
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


    //std::cout << "Copy from " << initial_word_id << " to " << final_word_id << std::endl;

    const uint8_t initial_read_offset = offset % 64;
    const uint8_t initial_read_length =
        std::min(bitcount, uint64_t(64) - initial_read_offset);

    rstr((*vector_)[initial_word_id]);
    writeInfix((*vector_)[initial_word_id],
               initial_read_offset,
               initial_read_length);

    if(final_word_id > initial_word_id) {
      reader_type reader(vector_begin_ + initial_word_id + 1, // skip first
                         vector_begin_ + final_word_id);      // skip last
      if(write_used > 0) {
        uint64_t read_word;
        while(!reader.empty()) {
          reader >> read_word;
          rstr(read_word);
          wstr((write_word | (read_word >> write_used)));
          writer_ << (write_word | (read_word >> write_used));
          write_word = (read_word << write_not_used);
        }
      } else {
        for(;!reader.empty(); ++reader) {rstr(*reader); writer_ << *reader; wstr(*reader);}
      }

      constexpr uint8_t final_read_offset = 0;
      const uint8_t final_read_length =
          (bitcount - initial_read_length + 63) % 64 + 1;
//          bitcount - initial_read_length - 64 * (final_word_id - initial_word_id - 1);

      rstr((*vector_)[final_word_id]);
      writeInfix((*vector_)[final_word_id],
                 final_read_offset,
                 final_read_length);
    }

    //std::cout << "CURRENT " << std::endl;
    rstr(write_word);
  }

  inline void finish() {
    writer_.finish();
  }


};

/******************************************************************************/
