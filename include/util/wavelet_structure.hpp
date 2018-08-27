/*******************************************************************************
 * include/util/common.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <memory>
#include <vector>
#include <type_traits>
#include <typeinfo>

#include "arrays/bit_vectors.hpp"
#include "huffman/huff_codes.hpp"

class wavelet_structure;

class base_wavelet_structure {
    friend class wavelet_structure;
public:
    // Prevent accidental copies
    base_wavelet_structure(base_wavelet_structure const&) = delete;
    base_wavelet_structure& operator =(base_wavelet_structure const&) = delete;

    // Allow moving
    base_wavelet_structure(base_wavelet_structure&& other) = default;
    base_wavelet_structure& operator =(base_wavelet_structure&& other) = default;

    virtual ~base_wavelet_structure() = default;
private:
    base_bit_vectors bvs_;
    size_t text_size_;
    bool is_tree_;
    bool is_huffman_shaped_;
protected:
    base_wavelet_structure(base_bit_vectors&& bvs, bool is_tree, bool is_huffman_shaped)
    : bvs_(std::move(bvs))
    , is_tree_(is_tree)
    , is_huffman_shaped_(is_huffman_shaped)
    {
        if (bvs_.levels() > 0) {
            text_size_ = bvs_.level_bit_size(0);
        } else {
            text_size_ = 0;
        }
    }

    inline virtual std::vector<uint64_t> const& zeros() const {
        static std::vector<uint64_t> empty;
        return empty;
    }
    inline virtual void const* codes(std::type_info const&) const {
        assert(false);
        return nullptr;
    }
};
class wavelet_structure_tree: public base_wavelet_structure {
public:
    wavelet_structure_tree()
    : wavelet_structure_tree(base_bit_vectors()) {}
    wavelet_structure_tree(base_bit_vectors&& bvs)
    : base_wavelet_structure(std::move(bvs), true, false) { }
};

class wavelet_structure_matrix: public base_wavelet_structure {
public:
    wavelet_structure_matrix()
    : wavelet_structure_matrix(base_bit_vectors(), {}) {}
    wavelet_structure_matrix(base_bit_vectors&& bvs,
                             std::vector<uint64_t>&& zeros)
    : base_wavelet_structure(std::move(bvs), false, false)
    , zeros_(std::move(zeros)) { }
private:
    std::vector<uint64_t> zeros_;
    inline virtual std::vector<uint64_t> const& zeros() const override {
        return zeros_;
    }
};

template<typename AlphabetType>
class wavelet_structure_tree_huffman: public base_wavelet_structure {
public:
    wavelet_structure_tree_huffman(canonical_huff_codes<AlphabetType, true>&& codes)
    : wavelet_structure_tree_huffman(base_bit_vectors(), std::move(codes)) { }
    wavelet_structure_tree_huffman(base_bit_vectors&& bvs,
                                   canonical_huff_codes<AlphabetType, true>&& codes)
    : base_wavelet_structure(std::move(bvs), true, true)
    , codes_(std::move(codes)) { }
private:
    canonical_huff_codes<AlphabetType, true> codes_;
    inline virtual void const* codes([[maybe_unused]] std::type_info const& type_info) const override {
        assert(type_info == typeid(canonical_huff_codes<AlphabetType, true>));
        return &codes_;
    }
};

template<typename AlphabetType>
class wavelet_structure_matrix_huffman: public base_wavelet_structure {
public:
    wavelet_structure_matrix_huffman(canonical_huff_codes<AlphabetType, false>&& codes)
    : wavelet_structure_matrix_huffman(base_bit_vectors(), {}, std::move(codes)) { }
    wavelet_structure_matrix_huffman(base_bit_vectors&& bvs,
                                     std::vector<uint64_t>&& zeros,
                                     canonical_huff_codes<AlphabetType, false>&& codes)
    : base_wavelet_structure(std::move(bvs), false, true)
    , codes_(std::move(codes))
    , zeros_(std::move(zeros)) { }
private:
    canonical_huff_codes<AlphabetType, false> codes_;
    std::vector<uint64_t> zeros_;

    inline virtual std::vector<uint64_t> const& zeros() const override {
        return zeros_;
    }

    inline virtual void const* codes([[maybe_unused]] std::type_info const& type_info) const override {
        assert(type_info == typeid(canonical_huff_codes<AlphabetType, false>));
        return &codes_;
    }
};

class wavelet_structure {
  wavelet_structure() = default;
public:
  template<typename concrete_wavelet_structure>
  wavelet_structure(concrete_wavelet_structure&& structure) {
    data_ = std::make_unique<concrete_wavelet_structure>(std::move(structure));
  }

  // Prevent accidental copies
  wavelet_structure(wavelet_structure const&) = delete;
  wavelet_structure& operator =(wavelet_structure const&) = delete;

  // Allow moving
  wavelet_structure(wavelet_structure&& other) = default;
  wavelet_structure& operator =(wavelet_structure&& other) = default;

  inline uint64_t levels() const {
    assert(bool(data_));
    return data_->bvs_.levels();
  }

  inline base_bit_vectors const& bvs() const {
    assert(bool(data_));
    return data_->bvs_;
  }

  inline std::vector<uint64_t> const& zeros() const {
    assert(bool(data_));
    return data_->zeros();
  }

  template<typename AlphabetType, bool is_tree>
  inline canonical_huff_codes<AlphabetType, is_tree> const& codes() const {
    using Ty = canonical_huff_codes<AlphabetType, is_tree>;
    assert(bool(data_));
    auto ptr = data_->codes(typeid(Ty));
    return *((Ty const*) ptr);
  }

  inline bool is_tree() const {
    assert(bool(data_));
    return data_->is_tree_;
  }

  inline bool is_huffman_shaped() const {
    assert(bool(data_));
    return data_->is_huffman_shaped_;
  }

  inline size_t text_size() const {
    assert(bool(data_));
    return data_->text_size_;
  }
private:
  std::unique_ptr<base_wavelet_structure> data_;
}; // class wavelet_structure

/******************************************************************************/
