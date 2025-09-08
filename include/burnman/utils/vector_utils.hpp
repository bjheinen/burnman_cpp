/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_VECTOR_UTILS_INCLUDED
#define BURNMAN_UTILS_VECTOR_UTILS_INCLUDED

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace utils {

  /**
  * @brief Helper to compute strides for mapping N-D grid to flat vector.
  *
  * For a grid with shape {2, 4}, the strides will be: {4, 1}.
  * The flat index can then be computed from the strides from a 2D index:
  *   {i0, i1} --> i0 * strides[0] + i1 * strides[1] --> i_flat
  * @see `flatten_index'
  * */
  inline std::vector<std::size_t> compute_strides(
    const std::vector<std::size_t>& shape
  ) {
    // Start with vector on ones (last stride always 1)
    std::vector<std::size_t> strides(shape.size(), 1);
    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(shape.size()) - 2; i >=0; --i) {
      std::size_t next_idx = static_cast<std::size_t>(i + 1);
      strides[static_cast<std::size_t>(i)] = strides[next_idx] * shape[next_idx];
    }
    return strides;
  }

  /**
  * @brief Maps N-D grid indices to flat vector index based on strides
  *
  * @see `compute_strides'
  */
  inline std::size_t flatten_index(
    const std::vector<std::size_t>& indices,
    const std::vector<std::size_t>& strides
  ) {
    if (indices.size() != strides.size()) {
      throw std::invalid_argument("Indices and strides must have the same size!");
    }
    std::size_t flat_index = 0;
    for (std::size_t i = 0; i < indices.size(); ++i) {
      flat_index += indices[i] * strides[i];
    }
    return flat_index;
  }

  /**
   * @brief Maps flat index back to N-D grid index given strides.
   *
   * @see `compute_strides'
   */
  inline std::vector<std::size_t> map_index(
    std::size_t flat_index,
    const std::vector<std::size_t>& strides
  ) {
    std::vector<std::size_t> indices(strides.size());
    for (std::size_t dim = 0; dim < strides.size(); ++dim) {
      indices[dim] = flat_index / strides[dim];
      flat_index %= strides[dim];
    }
  }

} // namespace utils

#endif // BURNMAN_UTILS_VECTOR_UTILS_INCLUDED
