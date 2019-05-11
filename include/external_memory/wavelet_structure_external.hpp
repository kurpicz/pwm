/*******************************************************************************
 * include/util/common.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>
#include "stxxl_helper.hpp"
#include "construction/wavelet_structure.hpp"
#include "util/debug_assert.hpp"

class wavelet_structure_external_writer;
class wavelet_structure_external_factory;

class wavelet_structure_external {
public:
  friend class wavelet_structure_external_writer;
  friend class wavelet_structure_external_factory;

  using bv_type = stxxlvector<uint64_t>;
  using zeros_type = std::vector<uint64_t>;
  using hists_type = std::vector<std::vector<uint64_t>>;

private:
  uint64_t text_size_;
  uint64_t levels_;

  bool is_tree_;
  bool save_zeros_;
  bool save_histograms_;
  bool is_huffman_shaped_;

//  bv_type metadata_;
  bv_type bvs_;

  uint64_t data_size_;
  std::vector<uint64_t> level_sizes_;
  std::vector<uint64_t> level_offsets_;

  zeros_type zeros_;
  hists_type histograms_;

  wavelet_structure_external(uint64_t text_size,
                             uint64_t levels,
                             bool is_tree,
                             bool save_zeros,
                             bool save_histograms,
                             bool is_huffman_shaped,
                             std::string em_name,
                             unsigned em_dir)
      : text_size_(text_size),
        levels_(levels),
        is_tree_(is_tree),
        save_zeros_(save_zeros),
        save_histograms_(save_histograms),
        is_huffman_shaped_(is_huffman_shaped),
//        metadata_(stxxl_files::getVectorPermanent<bv_type>(
//            em_dir, em_name + ".meta")),
        bvs_(stxxl_files::getVectorPermanent<bv_type>(
            em_dir, em_name + ".bvs")) {

    if (is_huffman_shaped_) std::abort(); //TODO: implement

    // make sure disk files are clean
    bvs_.clear();
//    metadata_.clear();

    data_size_ = 0;
    level_sizes_.resize(levels_, text_size_);
    level_offsets_.reserve(levels_ + 1);
    for (uint64_t level = 0; level < levels_; ++level) {
      level_offsets_.push_back(data_size_);
      data_size_ += (level_sizes_[level] + 63) / 64;
    }
    level_offsets_.push_back(data_size_);
    bvs_.resize(data_size_);

    if (save_zeros_) {
      zeros_.resize(levels, 0);
    } else {
      assert(is_tree_);
    }

    if (save_histograms_){
      uint64_t histogram_length = 1;
      histograms_.resize(levels + 1);
      for (uint64_t i = 0; i < levels + 1; i++) {
        histograms_[i].resize(histogram_length, 0);
        histogram_length *= 2;
      }
    }
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
    return level_sizes_;
  }

  inline std::vector<uint64_t> const& level_offsets() const {
    return level_offsets_;
  }

  inline std::vector<uint64_t> const& zeros() const {
    static std::vector<uint64_t> empty;
    if(!save_zeros_) return empty;
    return zeros_;
  }

  inline std::vector<std::vector<uint64_t>> const& histograms() const {
    static std::vector<std::vector<uint64_t>> empty;
    if(!save_histograms_) return empty;
    return histograms_;
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

  wavelet_structure getInternalStructure() {
    if(is_huffman_shaped_) std::abort();

    if(is_tree_) {
      return wavelet_structure_tree(getInternalBitvectors());
    } else {
      return wavelet_structure_matrix(getInternalBitvectors(),
                                      getInternalZeros());
    }
  }
};

class wavelet_structure_external_writer {
public:
  inline static auto& zeros(wavelet_structure_external& w) {
    return w.zeros_;
  }

  inline static auto& histograms(wavelet_structure_external& w) {
    return w.histograms_;
  }

  inline static auto& bvs(wavelet_structure_external& w) {
    return w.bvs_;
  }
};


class wavelet_structure_external_factory {
private:
  bool is_tree;
  bool save_zeros;
  bool save_hists;
  bool is_huff;
public:
  using self_type = wavelet_structure_external_factory;

  wavelet_structure_external_factory(bool tree)
      : is_tree(tree),
        save_zeros(!tree),
        save_hists(false),
        is_huff(false) {}

  wavelet_structure_external construct(uint64_t text_size,
                                       uint64_t levels,
                                       std::string em_name,
                                       unsigned em_dir) {

    return wavelet_structure_external(text_size,
                                      levels,
                                      is_tree,
                                      save_zeros,
                                      save_hists,
                                      is_huff,
                                      em_name,
                                      em_dir);
  }

  self_type zeros(bool value = true) {
    save_zeros = value || !is_tree;
    return *this;
  }

  self_type noZeros(bool value = true) {
    return zeros(!value);
  }


  self_type histograms(bool value = true) {
    save_hists = value;
    return *this;
  }

  self_type noHistograms(bool value = true) {
    return histograms(!value);
  }


  self_type huffman(bool value = true) {
    is_huff = value;
    if(is_huff) std::abort(); //TODO: implement
    return *this;
  }

  self_type noHuffman(bool value = true) {
    return huffman(!value);
  }
};

/******************************************************************************/
