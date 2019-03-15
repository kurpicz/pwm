#pragma once

template <typename char_t>
struct ps_ip_sort {
  const uint64_t n_;
  const uint64_t b_;

  struct mem_block {
    char_t *access_;
    uint64_t size_;

    inline char_t &operator[](const uint64_t idx) const {
      return access_[idx];
    }

    inline void push(const char_t c) {
      access_[size_++] = c;
    }
  };

  using q_type = std::queue<mem_block>;

  char_t *buffA_;
  char_t *buffB_;
  char_t *buffC_;
  char_t *buffD_;

  q_type q0_, q1_;
  std::vector<q_type *> qs_;

  std::queue<mem_block> q_free_;

  ps_ip_sort(char_t *text, const uint64_t n)
      : n_(n), b_(std::ceil(std::sqrt(n_))) {

    qs_.push_back(&q0_);
    qs_.push_back(&q1_);

    buffA_ = static_cast<char_t *>(malloc(b_ * sizeof(char_t)));
    buffB_ = static_cast<char_t *>(malloc(b_ * sizeof(char_t)));
    buffC_ = static_cast<char_t *>(malloc(b_ * sizeof(char_t)));

    q_free_.push(mem_block {buffA_, 0});
    q_free_.push(mem_block {buffB_, 0});
    q_free_.push(mem_block {buffC_, 0});

    const uint64_t last_block_size = n_ % b_;
    const uint64_t n_skip_ = n_ - last_block_size;

    for (uint64_t i = 0; i < n_skip_; i += b_) {
      q0_.push(mem_block {&text[i], b_} );
    }

    if (last_block_size > 0) {
      buffD_ = static_cast<char_t *>(malloc(b_ * sizeof(char_t)));
      for (uint64_t i = 0; i < last_block_size; ++i) {
        buffD_[i] = text[n_skip_ + i];
      }
      q0_.push(mem_block { buffD_, last_block_size });
    }
    else buffD_ = nullptr;
  }

  inline void sort(const uint8_t bit) {
    char_t mask = char_t(1ULL << (sizeof(char_t) * 8 - bit - 1));
    const uint64_t zero_blocks = q0_.size();
    const uint64_t one_blocks = q1_.size();

    mem_block b0_ = q_free_.front();
    q_free_.pop();
    mem_block b1_ = q_free_.front();
    q_free_.pop();

    for (uint64_t i = 0; i < zero_blocks; ++i) {
      auto block = q0_.front();
      auto blocksize = block.size_;
      block.size_ = 0;
      q0_.pop();
      q_free_.push(block);
      for (uint64_t j = 0; j < blocksize; ++j) {
        char_t character = block[j];
        std::cout << std::bitset<8>(character & mask).to_string() << std::endl;
        if (character & mask) {
          b1_.push(character);
          if (b1_.size_ == b_) {
            q1_.push(b1_);
            b1_ = q_free_.front();
            q_free_.pop();
          }
        }
        else {
          b0_.push(character);
          if (b0_.size_ == b_) {
            q0_.push(b0_);
            b0_ = q_free_.front();
            q_free_.pop();
          }
        }
      }
    }

    for (uint64_t i = 0; i < one_blocks; ++i) {
      auto block = q1_.front();
      auto blocksize = block.size_;
      block.size_ = 0;
      q1_.pop();
      q_free_.push(block);
      for (uint64_t j = 0; j < blocksize; ++j) {
        char_t character = block[j];
        if (character & mask) {
          b1_.push(character);
          if (b1_.size_ == b_) {
            q1_.push(b1_);
            b1_ = q_free_.front();
            q_free_.pop();
          }
        }
        else {
          b0_.push(character);
          if (b0_.size_ == b_) {
            q0_.push(b0_);
            b0_ = q_free_.front();
            q_free_.pop();
          }
        }
      }
    }

    q0_.push(b0_);
    q1_.push(b1_);
  }

  inline void print() {
    std::cout << "START PRINT n = " << n_ << ", b = " << b_ << std::endl;
    for (auto q_ : qs_) {
      for (uint64_t i = 0; i < q_->size(); ++i) {
        mem_block top = q_->front();
        q_->pop();

        for (uint64_t j = 0; j < top.size_; ++j) {
          auto character = top.access_[j];
          std::cout << std::bitset<8>(character).to_string();
          if (32 < character && character < 127) std::cout << "  ---  " << character;
          std::cout << std::endl;
        }
        q_->push(top);
      }
    }
  }

  ~ps_ip_sort() {
    delete buffA_;
    delete buffB_;
    delete buffC_;
    delete buffD_;
  }

};
