/*******************************************************************************
 * include/construction/pc_dd_fe/external_merge.hpp
 *
 * Copyright (C) 2019 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/


#pragma once

#include <deque>
#include "util/print.hpp"

#ifdef EXT_DD_CTX_VERBOSE
#define EXT_MERGE_VERBOSE EXT_DD_CTX_VERBOSE
#else
#define EXT_MERGE_VERBOSE if constexpr (false) atomic_out()
#endif

class em_merge_allocator {
private:
  constexpr static uint64_t default_blocks = 1024;
  uint64_t * data_;
  std::deque<uint64_t *> free_;
  uint64_t words_per_block_;

public:
  em_merge_allocator(const uint64_t min_blocks, uint64_t max_bytes) {
    const uint64_t max_words = std::max(min_blocks, max_bytes / 8);
    const uint64_t blocks = std::min(default_blocks, max_words);
    words_per_block_ = max_words / blocks;

    data_ = (uint64_t *) malloc(8ULL * words_per_block_ * blocks);
    for (uint64_t b = 0; b < blocks; ++b) {
      free_.push_back(&data_[b * words_per_block_]);
    }
  }

  ~em_merge_allocator() {
    delete data_;
  }

  inline uint64_t * get_block() {
    uint64_t * result = free_.front();
    free_.pop_front();
    return result;
  }

  inline void free_block(uint64_t * block) {
    free_.push_front(block);
  }

  inline uint64_t words_per_block() const {
    return words_per_block_;
  }

  inline bool full() const {
    return free_.empty();
  }

  inline uint64_t remaining_capacity_in_words() const {
    return free_.size() * words_per_block_;
  }

  em_merge_allocator(const em_merge_allocator&) = delete;
  em_merge_allocator& operator=(const em_merge_allocator&) = delete;
};


class em_merge_stack {
private:
  const uint64_t words_per_block_;
  em_merge_allocator& allocator_;
  std::deque<uint64_t *> blocks_;

  uint64_t * top_;
  uint64_t top_idx_;

  uint64_t * bot_;
  uint64_t bot_idx_;

  uint64_t bits_;
  std::deque<uint64_t> hist_;

  inline static uint64_t prefix(const uint64_t word, const uint64_t len) {
    const uint64_t mask = (0ULL - 1) << (64 - len);
    return word & mask;
  }

  inline static uint64_t suffix(const uint64_t word, const uint64_t len) {
    const uint64_t mask = (0ULL - 1) >> (64 - len);
    return word & mask;
  }

  inline static auto bitstring(const uint64_t word, std::string comment = "") {
    return std::bitset<64>(word).to_string() + " " + comment;
  }

public:
  em_merge_stack(em_merge_allocator &allocator)
      : words_per_block_(allocator.words_per_block()),
        allocator_(allocator),
        top_(nullptr), top_idx_(words_per_block_),
        bot_(nullptr), bot_idx_(0),
        bits_(0) {}

  inline uint64_t size_in_bits() const {
    return bits_;
  }

  template <typename reader_type>
  inline void to_memory(reader_type &reader, const uint64_t bits) {
    bits_ += bits;
    hist_.push_back(bits);

    const auto words = div64(bits + 63);
    for (uint64_t i = 0; i < words; ++i) {
      push(*reader);
      ++reader;
    }
  }

  template <typename writer_type>
  inline void to_em(writer_type &writer, uint64_t shift) {
    uint64_t initial_word = (shift == 0) ? 0ULL : prefix(*writer, shift);
    while(!hist_.empty()) {
      const uint64_t next_bits = hist_.front();
      const uint64_t words = div64(next_bits + 63);
      hist_.pop_front();

      EXT_MERGE_VERBOSE << "BITS " << next_bits << " WORDS " << words << " SHIFT " << shift << "\n";

      if (PWM_UNLIKELY(shift == 0)) {
        // all full words are perfectly aligned with the em vector
        // write full words
        for (uint64_t i = 0; i < (words - 1); ++i) {
          const uint64_t word = pop();
          EXT_MERGE_VERBOSE << bitstring(word, "POP -> WRITE (FULL ALIGN)\n");
          writer << word;
        }

        shift = mod64(next_bits);
        if (PWM_UNLIKELY(shift == 0)) {
          // last word is also full word
          const uint64_t word = pop();
          EXT_MERGE_VERBOSE << bitstring(word, "POP -> WRITE (FULL ALIGN AFTER FULL ALIGN)\n");
          writer << word;
          initial_word = 0ULL;
        }
        else {
          // last word is prefix
          initial_word = pop();
          EXT_MERGE_VERBOSE << bitstring(initial_word, "POP -> INITIAL WORD AFTER ALIGN\n");
        }
      }
      else {
        const uint64_t inverse_shift = 64 - shift;
        for (uint64_t i = 0; i < (words - 1); ++i) {
          const uint64_t word = pop();
          writer << (initial_word | (word >> shift));
          EXT_MERGE_VERBOSE << bitstring((initial_word | (word >> shift)), "POP -> WRITE (NO ALIGN)\n");
          initial_word = (word << inverse_shift);
          EXT_MERGE_VERBOSE << bitstring(initial_word, "    -> INITIAL WORD AFTER NO ALIGN\n");
        }

        const uint64_t last_bits = mod64(next_bits + 63) + 1;
        if (last_bits >= inverse_shift) {
          const uint64_t word = pop();
          writer << (initial_word | (word >> shift));
          EXT_MERGE_VERBOSE << bitstring((initial_word | (word >> shift)), "POP -> WRITE (NO ALIGN)\n");
          initial_word = (word << inverse_shift);
          shift = last_bits - inverse_shift;
          EXT_MERGE_VERBOSE << ((shift > 0) ? bitstring(initial_word, "    -> INITIAL WORD AFTER NO ALIGN\n") : "");
        }
        else {
          const uint64_t word = pop();
          initial_word = (initial_word | (word >> shift));
          EXT_MERGE_VERBOSE << bitstring(initial_word, "POP -> INITIAL WORD AFTER NO ALIGN\n");
          shift = mod64(shift + next_bits);
        }
      }
    }

    if (PWM_LIKELY(shift > 0)) {
      EXT_MERGE_VERBOSE << bitstring(initial_word) << "\n";
      writer << (initial_word | suffix(*writer, 64 - shift));
    }
    clear();
  }

private:
  // push back
  inline void push(const uint64_t value) {
    if (PWM_UNLIKELY(top_idx_ == words_per_block_)) {
      blocks_.push_back(allocator_.get_block());
      top_ = blocks_.back();
      top_[0] = value;
      top_idx_ = 1;
      if (PWM_UNLIKELY(blocks_.size() == 1)) {
        bot_ = top_;
        bot_idx_ = 0;
      }
    }
    else top_[top_idx_++] = value;
  }

  // pop front
  inline uint64_t pop() {
    if (PWM_UNLIKELY(bot_idx_ == words_per_block_)) {
      allocator_.free_block(bot_);
      blocks_.pop_front();
      bot_ = blocks_.front();
      bot_idx_ = 1;
      return bot_[0];
    }
    else return bot_[bot_idx_++];
  }

  inline void clear() {
    while (!blocks_.empty()) {
      allocator_.free_block(blocks_.back());
      blocks_.pop_back();
    }
    top_ = nullptr;
    top_idx_ = words_per_block_;
    bot_ = nullptr;
    bot_idx_ = 0;
    bits_ = 0;
    hist_.clear();
  }
};


class em_merger {
private:
  const uint64_t levels_;
  const uint64_t intervals_;
  em_merge_allocator allocator_;
  std::vector<em_merge_stack> buffers_;

public:
  em_merger(const uint64_t levels, const uint64_t bytes_memory)
      : levels_(levels),
        intervals_((1ULL << levels) - 1),
        allocator_(intervals_, bytes_memory) {
    buffers_.reserve(intervals_);
    for (uint64_t i = 0; i < intervals_; ++i) {
      buffers_.emplace_back(allocator_);
    }
  }

  inline auto &interval(const uint64_t idx) {
    return buffers_[idx];
  }

  // TODO: this could be done using a priority queue
  inline uint64_t largest_interval() const {
    uint64_t result = 0;
    uint64_t max_size = buffers_[0].size_in_bits();
    for (uint64_t i = 1; i < intervals_; ++i) {
      const uint64_t cur_size = buffers_[i].size_in_bits();
      result = (cur_size > max_size) ? i : result;
      max_size = std::max(max_size, cur_size);
    }
    return result;
  }

  inline uint64_t remaining_capacity_in_words() const {
    return allocator_.remaining_capacity_in_words();
  }
};


















