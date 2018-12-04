/*******************************************************************************
 * include/util/common.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>
#include "arrays/stxxl_helper.hpp"
#include "construction/wavelet_structure.hpp"
#include "util/debug_assert.hpp"

class wavelet_structure_external_writer;

class wavelet_structure_external {
public:
  friend class wavelet_structure_external_writer;
private:
  uint64_t levels_;
  uint64_t text_size_;
  std::vector<uint64_t> * level_sizes_;
  std::vector<uint64_t> * zeros_;
  std::vector<std::vector<uint64_t>> * histograms_;

  bool is_tree_;
  bool is_huffman_shaped_;
  bool save_zeros_;
  bool save_histograms_;

  stxxlvector<uint64_t> metadata_;
  stxxlvector<uint64_t> bvs_;
  std::vector<uint64_t> level_offsets_;
  uint64_t data_size_;

  wavelet_structure_external(std::vector<uint64_t> * level_sizes,
                             std::vector<uint64_t> * zeros,
                             std::vector<std::vector<uint64_t>> * histograms,
                             bool is_tree,
                             bool is_huffman_shaped,
                             unsigned em_dir,
                             std::string em_name)
      : levels_((level_sizes == nullptr) ? 0 : level_sizes->size()),
        text_size_((levels_ > 0) ? (*level_sizes)[0] : 0),
        level_sizes_((level_sizes == nullptr) ? new std::vector<uint64_t>() : level_sizes),
        zeros_(zeros),
        histograms_(histograms),
        is_tree_(is_tree),
        is_huffman_shaped_(is_huffman_shaped),
        save_zeros_(zeros_ != nullptr),
        save_histograms_(histograms_ != nullptr) {
    if(save_zeros_) DCHECK(zeros_->size() == levels_);
    if(save_histograms_) {
      DCHECK(histograms_->size() == levels_ + 1);
      (*histograms_)[0][0] = text_size_;
    }

    metadata_ = stxxl_files::getVectorPermanent<stxxlvector<uint64_t>>(em_dir, em_name + ".meta");
    bvs_ = stxxl_files::getVectorPermanent<stxxlvector<uint64_t>>(em_dir, em_name + ".bvs");

    bvs_.clear();
    metadata_.clear();

    level_offsets_.reserve(levels_ + 1);
    data_size_ = 0;
    for (uint64_t level = 0; level < levels_; ++level) {
      level_offsets_.push_back(data_size_);
      data_size_ += ((*level_sizes_)[level] + 63) / 64;
    }
    level_offsets_.push_back(data_size_);
    bvs_.resize(data_size_);
  }

  bit_vectors<> getInternalBitvectors() {
    if(is_huffman_shaped_) std::abort();

    bit_vectors<> result(levels_, text_size_);
    for(uint64_t level = 0; level < levels_; level++) {
      auto intLevel = result[level];
      auto extLevel = (*this)[level];
      for(uint64_t datapos = 0; datapos < intLevel.size(); datapos++) {
        intLevel[datapos] = extLevel[datapos];
      }
    }
    return result;
  }

  std::vector<uint64_t> getInternalZeros() {
    if(is_huffman_shaped_) std::abort();

    std::vector<uint64_t> result(levels_);
    for(uint64_t level = 0; level < levels_; level++) {
      result[level] = zeros()[level];
    }
    return result;
  }


public:
  // Prevent accidental copies
  wavelet_structure_external(wavelet_structure_external const&) = delete;
  wavelet_structure_external& operator=(wavelet_structure_external const&) = delete;

  // Allow moving
  wavelet_structure_external(wavelet_structure_external&& other) = default;
  wavelet_structure_external& operator=(wavelet_structure_external&& other) = default;

  inline bool is_tree() const {
    return is_tree_;
  }

  inline bool has_zeros() const {
    return save_zeros_;
  }

  inline bool has_histograms() const {
    return save_histograms_;
  }

  inline uint64_t levels() const {
    return levels_;
  }

  inline uint64_t text_size() const {
    return text_size_;
  }

  inline std::vector<uint64_t> const& level_sizes() const {
    return *level_sizes_;
  }

  inline std::vector<uint64_t> const& level_offsets() const {
    return level_offsets_;
  }

  inline std::vector<uint64_t> const& zeros() const {
    static std::vector<uint64_t> empty;
    if(!save_zeros_) return empty;
    return *zeros_;
  }

  inline std::vector<std::vector<uint64_t>> const& histograms() const {
    static std::vector<std::vector<uint64_t>> empty;
    if(!save_histograms_) return empty;
    return *histograms_;
  }

  inline const stxxlvector_offset<uint64_t>
  operator[](const uint64_t level) const {
    return stxxlvector_offset<uint64_t>(bvs_, level_offsets_[level]);
  }

  inline stxxlvector_offset<uint64_t> operator[](const uint64_t level) {
    return stxxlvector_offset<uint64_t>(bvs_, level_offsets_[level]);
  }

  void saveMetaData() {
    // TODO: implement
  }

  ~wavelet_structure_external() {
    delete level_sizes_;
    delete zeros_;
    delete histograms_;
  }

  wavelet_structure getInternalStructure() {
    if(is_huffman_shaped_) std::abort();

    if(is_tree_) {
      return wavelet_structure_tree(getInternalBitvectors());
    } else {
      return wavelet_structure_matrix(getInternalBitvectors(),
                                      getInternalZeros());
    }
  }

  static wavelet_structure_external getNoHuffman(
      uint64_t size,
      uint64_t levels,
      bool is_tree,
      unsigned em_dir,
      std::string em_name) {
    std::vector<uint64_t> * level_sizes =
        new std::vector<uint64_t>(levels, size);
    std::vector<uint64_t> * zeros = nullptr;
    std::vector<std::vector<uint64_t>> * histograms = nullptr;

    // TODO: do not create histograms if not needed
    {
      uint64_t histogram_length = 1;
      histograms = new std::vector <std::vector<uint64_t>>(levels + 1);
      for (uint64_t i = 0; i < levels + 1; i++) {
        (*histograms)[i].resize(histogram_length, 0);
        histogram_length *= 2;
      }
    }
    if(!is_tree) {
      zeros = new std::vector<uint64_t>(levels, 0);
    }

    return wavelet_structure_external(level_sizes,
                                      zeros,
                                      histograms,
                                      is_tree,
                                      false,
                                      em_dir,
                                      em_name);
  }
};

class wavelet_structure_external_writer {
public:
  inline static auto& zeros(wavelet_structure_external& w) {
    return *(w.zeros_);
  }

  inline static auto& histograms(wavelet_structure_external& w) {
    return *(w.histograms_);
  }

  inline static auto& bvs(wavelet_structure_external& w) {
    return w.bvs_;
  }
};

/******************************************************************************/
