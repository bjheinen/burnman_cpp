/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_TYPES_HPP_INCLUDED
#define BURNMAN_UTILS_TYPES_HPP_INCLUDED

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <variant>

/**
 * Custom type alias for chemical formulae.
 */
using FormulaMap = std::unordered_map<std::string, double>;

// Overloaded operators to add and weight formulae
inline FormulaMap operator+(const FormulaMap& a, const FormulaMap& b) {
  FormulaMap result = a;
  for (const auto& [elem, count] : b) {
    result[elem] += count;
  }
  return result;
}

inline FormulaMap& operator+=(FormulaMap& a, const FormulaMap& b) {
  for (const auto& [elem, count] : b) {
    a[elem] += count;
  }
  return a;
}

inline FormulaMap operator*(const FormulaMap& a, double scalar) {
  FormulaMap result;
  for (const auto& [elem, count] : a) {
    result[elem] = scalar * count;
  }
  return result;
}

inline FormulaMap operator*(double scalar, const FormulaMap& a) {
  return a * scalar;
}

inline FormulaMap& operator*=(FormulaMap& a, double scalar) {
  for (auto& [_, count] : a) {
    count *= scalar;
  }
  return a;
}

/**
 * Enum used to define EOS Type
 */
enum class EOSType {
  Auto, // Used to set EOS from params
  Custom, // Used when user passes custom EOSType
  Vinet,
  BM3,
  BM2,
  MGD2,
  MGD3,
  SLB2,
  SLB3,
  SLB3Conductive
};

/**
 * Enum used to define averaging scheme type
 */
enum class AveragingType {
  Voigt,
  Reuss,
  VRH,
  HashinShtrikmanLower,
  HashinShtrikmanUpper,
  HashinShtrikman
};

/**
 * Enum used to define fraction type
 */
enum class FractionType {
  Molar,
  Mass,
  Volume
};

/**
 * Forward declarations
 */
class Mineral;

/**
 * Type aliases for constructing solid solution models.
 */
using MineralFormulaPair = std::pair<Mineral, std::string>;
using PairedEndmemberList = std::vector<MineralFormulaPair>;

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
    : shape(shape), strides(utils::compute_strides(shape), size(compute_size(shape)))
  {
    this->data.resize(this->size);
  }

  /**
   * @brief Convenience function for data access
   */
  std::vector<T>& data() {
    return this->data;
  }

  /**
   * @brief Convenience function for read-only data access.
   */
  const std::vector<T>& data() const {
    return this->data;
  }

  /**
   * @brief Returns NDArray shape.
   */
  const std::vector<std::size_t>& shape() const {
    return this->shape;
  }

  /**
   * @brief Returns data size (total number of elements).
   */
  const std::size_t size() const {
    return this->data.size();
  }

  /**
   * @brief Sets data shape.
   *
   * @note Clears any existing data. To be used after default constructing.
   */
  void set_shape(const std::vector<std::size_t>& shape) {
    this->shape = shape;
    this->strides = utils::compute_strides(this->shape);
    this->size = compute_size(this->shape);
    // Clear data
    this->data.clear();
    this->data.resize(this->size);
  }

  // Indexing
  // Flat indexing
  // mutable, with []
  T& operator[](std::size_t flat_index) {
    return this->data[flat_index];
  }
  // const, with []
  const T& operator[](std::size_t flat_index) const {
    return this->data[flat_index];
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
    return this->data[utils::flatten_index(indices, this->strides)];
  }
  // const, with vector
  const T& operator()(const std::vector<std::size_t>& indices) const {
    return this->data[utils::flatten_index(indices, this->strides)];
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
  std::vector<T> data;
  std::vector<std::size_t> shape;
  std::vector<std::size_t> strides;
  std::size_t size;

  static std::size_t compute_size(const std::vector<std::size_t>& new_shape) {
    std::size_t new_size = 1;
    for (std::size_t dim : new_shape) {
      new_size *= dim;
    }
    return new_size;
  }

}

#endif // BURNMAN_UTILS_TYPES_HPP_INCLUDED
