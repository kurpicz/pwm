/*******************************************************************************
 * include/util/stats.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <chrono>
#include <ostream>
#include <sstream>

#include <stxxl/bits/io/iostats.h>

#ifdef MALLOC_COUNT
#include "benchmark/malloc_count.h"
#endif // MALLOC_COUNT

template <bool stxxl_stats>
struct statistics;

template <bool stxxl_stats>
struct stxxl_stats_type;

template <>
struct stxxl_stats_type<true> {
  friend class statistics<true>;

  std::string title_;

  stxxl::stats * get_stats_;
  stxxl::stats_data start_stats_;
  stxxl::stats_data total_stats_;

  stxxl_stats_type(std::string title)
      : title_(title),
        get_stats_(stxxl::stats::get_instance()) {}

  inline void start() {
    start_stats_ = stxxl::stats_data(*get_stats_);
  }

  inline void finish() {
    total_stats_ = stxxl::stats_data(*get_stats_) - start_stats_;
  }

  std::string to_string() const {
    std::stringstream result;
    const stxxl::stats_data &s = total_stats_;

    result <<        title_ << "r_bytes=" << s.get_read_volume();
    result << " " << title_ << "r_blocks=" << s.get_reads();
    result << " " << title_ << "w_bytes=" << s.get_written_volume();
    result << " " << title_ << "w_blocks=" << s.get_writes();

    result << " " << title_ << "r_bytes_cached=" << s.get_cached_read_volume();
    result << " " << title_ << "r_blocks_cached=" << s.get_cached_reads();
    result << " " << title_ << "w_bytes_cached=" << s.get_cached_written_volume();
    result << " " << title_ << "w_blocks_cached=" << s.get_cached_writes();

    result << " " << title_ << "r_time=" << uint64_t(s.get_read_time() * 1000);
    result << " " << title_ << "w_time=" << uint64_t(s.get_write_time() * 1000);
    result << " " << title_ << "r_wait_time="
                  << uint64_t(s.get_wait_read_time() * 1000);
    result << " " << title_ << "w_wait_time="
                  << uint64_t(s.get_wait_write_time() * 1000);
    result << " " << title_ << "rw_wait_time="
                  << uint64_t(s.get_io_wait_time() * 1000);
    result << " " << title_ << "r_parallel_time="
                  << uint64_t(s.get_pread_time() * 1000);
    result << " " << title_ << "w_parallel_time="
                  << uint64_t(s.get_pwrite_time() * 1000);
    result << " " << title_ << "rw_parallel_time="
                  << uint64_t(s.get_pio_time() * 1000);

    return result.str();
  }
};

template <>
struct stxxl_stats_type<false> {
  stxxl_stats_type(std::string) {}
  inline void start() {}
  inline void finish() {}
};

template <bool stxxl_stats>
class statistics {
private:

  const statistics *parent_;
  std::string title_;

  decltype(std::chrono::high_resolution_clock::now()) start_time_;
  uint64_t total_time_;

  int64_t start_memory_;
  int64_t total_memory_;

  stxxl_stats_type<stxxl_stats> stxxl_stats_;

  std::vector<statistics> phases_;
  bool phase_running_ = false;

  inline statistics(std::string title, const statistics *parent) :
      parent_(parent), title_(title), stxxl_stats_(title) {
    if (parent != nullptr)
      start_memory_ = parent->start_memory_;
  }

public:

  inline statistics() : statistics("", nullptr) {}

  inline void start() {
    stxxl_stats_.start();
    #ifdef MALLOC_COUNT
    malloc_count_reset_peak();
    if(parent_ == nullptr)
      start_memory_ = malloc_count_peak();
    #endif // MALLOC_COUNT
    start_time_ = std::chrono::high_resolution_clock::now();
  }

  inline void finish() {
    if (phase_running_)
      phases_[phases_.size() - 1].finish();
    phase_running_ = false;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = end_time - start_time_;
    total_time_ =
        std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
    #ifdef MALLOC_COUNT
    total_memory_ = malloc_count_peak() - start_memory_;
    for (auto phase : phases_)
      total_memory_ = std::max(total_memory_, phase.total_memory_);
    #endif // MALLOC_COUNT
    stxxl_stats_.finish();
  }

  inline void phase(std::string title = "") {
    if (phase_running_)
      phases_[phases_.size() - 1].finish();

    phases_.push_back(statistics((title.empty()
        ? "Phase " + std::to_string(phases_.size() + 1) : title) + "_", this));
    phases_[phases_.size() - 1].start();
    phase_running_ = true;
  }

  uint64_t get_total_memory() const {
    return total_memory_;
  }

  uint64_t get_total_time() const {
    return total_time_;
  }

  bool operator < (const statistics& other) const {
    return (start_time_ < other.start_time_);
  }

  std::string to_string() const {
    std::stringstream result;

    result << (title_.empty() ? "median_" : title_) << "time=" << total_time_;

    #ifdef MALLOC_COUNT
    result << " " << title_ << "additional_memory=" << total_memory_;
    #endif // MALLOC_COUNT

    if constexpr (stxxl_stats)
      result << " " << stxxl_stats_.to_string();

    if (!phases_.empty())
      result << " phases="
             << phases_[0].title_.substr(0, phases_[0].title_.size() - 1);

    for (uint64_t i = 1; i < phases_.size(); ++i)
      result << ","
             << phases_[i].title_.substr(0, phases_[i].title_.size() - 1);

    for(auto phase : phases_) {
      result << " " << phase.to_string();
    }

    return result.str();
  }

};

template <bool stxxl_stats>
std::ostream &operator<<(std::ostream &os, statistics<stxxl_stats> const &m) {
  return os << m.to_string();
}