/*******************************************************************************
 * include/util/merge.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "../../../../../../../usr/include/c++/7/bitset"
#include "../../arrays/stxxl_helper.hpp"

template <typename vector_type>
class external_merger {
  using value_type = typename vector_type::value_type;
  using writer_type = typename vector_type::bufwriter_type;
  using reader_type = typename vector_type::bufreader_type;
  using iter_type = typename vector_type ::const_iterator;

  writer_type writer_;
  const vector_type& vector_;
  const iter_type vector_begin_;
  const iter_type vector_end_;

  uint64_t write_word = 0;
  uint8_t write_used = 0;
  uint8_t write_not_used = 64;

  std::vector<uint64_t> left_masks_vec;
  uint64_t * left_mask;

public:
  external_merger(const vector_type& input, vector_type& output)
      : writer_(output),
        vector_(input),
        vector_begin_(vector_.begin()),
        vector_end_(vector_.end()){
    static_assert(sizeof(value_type) == 8);

    for(uint8_t width = 0; width < 64; width++) {
      left_masks_vec.push_back(((0x1ULL << width) - 1) << (64 - width));
    }
    left_masks_vec.push_back(~uint64_t(0));
    left_mask = left_masks_vec.data();
  }

  inline void finishLevel() {
    if(write_used > 0) {
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

    // calculate indices
    const uint64_t initial_word_id = offset / 64;
    const uint64_t final_word_id = (offset + bitcount - 1) / 64;
    const uint8_t initial_read_offset = offset % 64;
    const uint8_t initial_read_length =
        std::min(bitcount, uint64_t(64) - initial_read_offset);
    const uint8_t final_read_length =
        (bitcount - initial_read_length + 63) % 64 + 1;

    // write initial word
    writeInfix(vector_[initial_word_id],
               initial_read_offset,
               initial_read_length);

    if(final_word_id > initial_word_id) {
      // use bufreader for large ranges
      if(bitcount > 8ULL * 128 * 1024 * 1024) {
        reader_type reader(vector_begin_ + initial_word_id + 1, // skip first
                           vector_begin_ + final_word_id);      // skip last
        if(write_used > 0) {
          uint64_t read_word;
          while(!reader.empty()) {
            reader >> read_word;
            writer_ << (write_word | (read_word >> write_used));
            write_word = (read_word << write_not_used);
          }
        } else {
          for(;!reader.empty(); ++reader) writer_ << *reader;
        }
      }
      // use index access for small ranges
      else {
        uint64_t current_word_id = initial_word_id + 1; //skip first
        if(write_used > 0) {
          uint64_t read_word;
          while(current_word_id < final_word_id) {
            read_word = vector_[current_word_id++];
            writer_ << (write_word | (read_word >> write_used));
            write_word = (read_word << write_not_used);
          }
        } else {
          for(;current_word_id < final_word_id; ++current_word_id)
            writer_ << vector_[current_word_id];
        }
      }

      // write final word
      writeInfix(vector_[final_word_id], 0, final_read_length);
    }
  }


  inline void finish() {
    writer_.finish();
  }


};



template <typename vector_type>
class external_buffered_merger {
  using value_type = typename vector_type::value_type;
  using writer_type = typename vector_type::bufwriter_type;
  using reader_type = typename vector_type::bufreader_type;
  using iter_type = typename vector_type ::const_iterator;

  writer_type writer_;
  const vector_type& vector_;
  const iter_type vector_begin_;
  const iter_type vector_end_;

  uint64_t write_word = 0;
  uint8_t write_used = 0;
  uint8_t write_not_used = 64;

  std::vector<uint64_t> left_masks_vec;
  uint64_t * left_mask;

public:
  external_buffered_merger(const vector_type& input, vector_type& output)
      : writer_(output),
        vector_(input),
        vector_begin_(vector_.begin()),
        vector_end_(vector_.end()){
    static_assert(sizeof(value_type) == 8);

    for(uint8_t width = 0; width < 64; width++) {
      left_masks_vec.push_back(((0x1ULL << width) - 1) << (64 - width));
    }
    left_masks_vec.push_back(~uint64_t(0));
    left_mask = left_masks_vec.data();
  }

  inline void finishLevel() {
    if(write_used > 0) {
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

    // calculate indices
    const uint64_t initial_word_id = offset / 64;
    const uint64_t final_word_id = (offset + bitcount - 1) / 64;
    const uint8_t initial_read_offset = offset % 64;
    const uint8_t initial_read_length =
        std::min(bitcount, uint64_t(64) - initial_read_offset);
    const uint8_t final_read_length =
        (bitcount - initial_read_length + 63) % 64 + 1;

    // write initial word
    writeInfix(vector_[initial_word_id],
               initial_read_offset,
               initial_read_length);

    if(final_word_id > initial_word_id) {
      // use bufreader for large ranges
      if(bitcount > 8ULL * 128 * 1024 * 1024) {
        reader_type reader(vector_begin_ + initial_word_id + 1, // skip first
                           vector_begin_ + final_word_id);      // skip last
        if(write_used > 0) {
          uint64_t read_word;
          while(!reader.empty()) {
            reader >> read_word;
            writer_ << (write_word | (read_word >> write_used));
            write_word = (read_word << write_not_used);
          }
        } else {
          for(;!reader.empty(); ++reader) writer_ << *reader;
        }
      }
        // use index access for small ranges
      else {
        uint64_t current_word_id = initial_word_id + 1; //skip first
        if(write_used > 0) {
          uint64_t read_word;
          while(current_word_id < final_word_id) {
            read_word = vector_[current_word_id++];
            writer_ << (write_word | (read_word >> write_used));
            write_word = (read_word << write_not_used);
          }
        } else {
          for(;current_word_id < final_word_id; ++current_word_id)
            writer_ << vector_[current_word_id];
        }
      }

      // write final word
      writeInfix(vector_[final_word_id], 0, final_read_length);
    }
  }


  inline void finish() {
    writer_.finish();
  }
};


/******************************************************************************/
