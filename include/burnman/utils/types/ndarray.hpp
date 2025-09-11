/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_TYPES_NDARRAY_HPP_INCLUDED
#define BURNMAN_UTILS_TYPES_NDARRAY_HPP_INCLUDED

#include <cstddef>
#include <initializer_list>
#include <vector>
#include "burnman/utils/vector_utils.hpp"

/**
 * @brief Flat vector storage with convenience mapping to grid indices.
 *
 * Arrays can be default (empty) constructed or constructed with a shape vector.
 * When default constructing, use `arr.set_shape(shape)' before setting any data.
 *
 * Use `arr.data()' to access the data vector directly.
 * Use `arr.shape()' for shape.
 * Use `arr.size()' for size.
 *
 * Indexing:
 *   Indexing can be done with a flat index, i, or grid indices {i, j, k, ...}.
 *   Overloads are provided for () and [] for flat indices, and () for grid indices.
 *   Initializer lists can also be used:
 *    arr({2, 3})
 *
 */
template<typename T>
struct NDArray {
 public:
  /**
   * @brief Default construct and empty NDArray.
   */
  NDArray() = default;

  /**
   * @brief Construct an NDArray with a given shape.
   */
  explicit NDArray(const std::vector<std::size_t>& shape)
    : shape_(shape), strides_(utils::compute_strides(shape), size_(compute_size(shape)))
  {
    this->storage_.resize(this->size_);
  }

  /**
   * @brief Convenience function for data access
   */
  std::vector<T>& data() {
    return this->storage_;
  }

  /**
   * @brief Convenience function for read-only data access.
   */
  const std::vector<T>& data() const {
    return this->storage_;
  }

  /**
   * @brief Returns NDArray shape.
   */
  const std::vector<std::size_t>& shape() const {
    return this->shape_;
  }

  /**
   * @brief Returns data size (total number of elements).
   */
  std::size_t size() const {
    return this->size_;
  }

  /**
   * @brief Sets data shape.
   *
   * @note Clears any existing data. To be used after default constructing.
   */
  void set_shape(const std::vector<std::size_t>& shape) {
    this->shape_ = shape;
    this->strides_ = utils::compute_strides(this->shape_);
    this->size_ = compute_size(this->shape_);
    // Clear data
    this->storage_.clear();
    this->storage_.resize(this->size_);
  }

  // Indexing
  // Flat indexing
  // mutable, with []
  T& operator[](std::size_t flat_index) {
    return this->storage_[flat_index];
  }
  // const, with []
  const T& operator[](std::size_t flat_index) const {
    return this->storage_[flat_index];
  }
  // mutable, with ()
  T& operator()(std::size_t flat_index) {
    return (*this)[flat_index];
  }
  // const, with ()
  const T& operator()(std::size_t flat_index) const {
    return (*this)[flat_index];
  }
  // Grid indexing
  // mutable, with vector
  T& operator()(const std::vector<std::size_t>& indices) {
    return this->storage_[utils::flatten_index(indices, this->strides_)];
  }
  // const, with vector
  const T& operator()(const std::vector<std::size_t>& indices) const {
    return this->storage_[utils::flatten_index(indices, this->strides_)];
  }
  // mutable, with initialiser list
  T& operator()(std::initializer_list<std::size_t> indices) {
    return (*this)(std::vector<std::size_t>(indices));
  }
  // const, with initialiser list
  const T& operator()(std::initializer_list<std::size_t> indices) const {
    return (*this)(std::vector<std::size_t>(indices));
  }

 private:
  std::vector<T> storage_;
  std::vector<std::size_t> shape_;
  std::vector<std::size_t> strides_;
  std::size_t size_;

  static std::size_t compute_size(const std::vector<std::size_t>& shape) {
    std::size_t size = 1;
    for (std::size_t dim : shape) {
      size *= dim;
    }
    return size;
  }

};

#endif // BURNMAN_UTILS_TYPES_NDARRAY_HPP_INCLUDED
