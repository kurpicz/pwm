#pragma once

#include <bitset>
#include <util/macros.hpp>

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
  char_t *buffLastBlock_;

  q_type q0_, q1_;

  std::queue<mem_block> q_free_;

  ps_ip_sort(char_t *text, const uint64_t n, const uint64_t b)
      : n_(n), b_(b) {

    buffA_ = static_cast<char_t *>(malloc(b_ * sizeof(char_t)));
    buffB_ = static_cast<char_t *>(malloc(b_ * sizeof(char_t)));
    buffC_ = static_cast<char_t *>(malloc(b_ * sizeof(char_t)));
    buffD_ = static_cast<char_t *>(malloc(b_ * sizeof(char_t)));

    q_free_.push(mem_block {buffA_, 0});
    q_free_.push(mem_block {buffB_, 0});
    q_free_.push(mem_block {buffC_, 0});
    q_free_.push(mem_block {buffD_, 0});

    const uint64_t last_block_size = n_ % b_; //might be 0
    const uint64_t n_skip_ = n_ - last_block_size;

    for (uint64_t i = 0; i < n_skip_; i += b_) {
      q0_.push(mem_block {&text[i], b_} );
    }

    if (last_block_size > 0) {
      buffLastBlock_ = static_cast<char_t *>(malloc(b_ * sizeof(char_t)));
      for (uint64_t i = 0; i < last_block_size; ++i) {
        buffLastBlock_[i] = text[n_skip_ + i];
      }
      q0_.push(mem_block { buffLastBlock_, last_block_size });
    }
    else buffLastBlock_ = nullptr;

//    std::cout << "\n\nINIT:  " << q_free_.size() << " " << q0_.size() << " " << q1_.size() << std::endl;
  }

  inline uint64_t get_zeros() const { return (q0_.size() == 0 ? 0 : (q0_.size() - 1) * b_ + q0_.back().size_); }
  inline uint64_t get_ones() const { return (q1_.size() == 0 ? 0 : (q1_.size() - 1) * b_ + q1_.back().size_); }

  uint64_t zero_idx_pop;
  uint64_t one_idx_pop;
  mem_block zero_block_pop;
  mem_block one_block_pop;

  mem_block zero_block_push;
  mem_block one_block_push;

  inline void start_level() {
//    std::cout << "START A: " << q_free_.size() << " " << q0_.size() << " " << q1_.size() << std::endl;
    const auto zeros = get_zeros();
    if (zeros > 0) {
      zero_idx_pop = 0;
      zero_block_pop = q0_.front();
    }
    if (n_ - zeros > 0) {
      one_idx_pop = 0;
      one_block_pop = q1_.front();
    }

    zero_block_push = q_free_.front();
    zero_block_push.size_ = 0;
    q_free_.pop();
    one_block_push = q_free_.front();
    one_block_push.size_ = 0;
    q_free_.pop();
//    std::cout << "START B: " << q_free_.size() << " " << q0_.size() << " " << q1_.size() << std::endl;
  }

  inline void finish_level() {
//    std::cout << "STOP A:  " << q_free_.size() << " " << q0_.size() << " " << q1_.size() << std::endl;
    if (zero_block_push.size_ > 0)
      q0_.push(zero_block_push);
    else
      q_free_.push(zero_block_push);

    if (one_block_push.size_ > 0)
      q1_.push(one_block_push);
    else
      q_free_.push(one_block_push);

//    std::cout << "STOP B:  " << q_free_.size() << " " << q0_.size() << " " << q1_.size() << std::endl;
  }

  inline char_t pop0() {
    if (PWM_UNLIKELY(zero_idx_pop == zero_block_pop.size_ - 1)) {
      const auto result = zero_block_pop[zero_idx_pop];
      zero_idx_pop = 0;
      q_free_.push(zero_block_pop);
      q0_.pop();
      zero_block_pop = q0_.front();
      return result;
    }
    return zero_block_pop[zero_idx_pop++];
  }

  inline char_t pop1() {
    if (PWM_UNLIKELY(one_idx_pop == one_block_pop.size_ - 1)) {
      const auto result = one_block_pop[one_idx_pop];
      one_idx_pop = 0;
      q_free_.push(one_block_pop);
      q1_.pop();
      one_block_pop = q1_.front();
      return result;
    }
    return one_block_pop[one_idx_pop++];
  }

  inline void push0(const char_t character) {
    zero_block_push.push(character);
    if (PWM_UNLIKELY(zero_block_push.size_ == b_)) {
      q0_.push(zero_block_push);
      zero_block_push = q_free_.front();
      zero_block_push.size_ = 0;
      q_free_.pop();
    }
  }

  inline void push1(const char_t character) {
    one_block_push.push(character);
    if (PWM_UNLIKELY(one_block_push.size_ == b_)) {
      q1_.push(one_block_push);
      one_block_push = q_free_.front();
      one_block_push.size_ = 0;
      q_free_.pop();
    }
  }

  ~ps_ip_sort() {
    delete buffA_;
    delete buffB_;
    delete buffC_;
    delete buffD_;
    delete buffLastBlock_;
  }
};
