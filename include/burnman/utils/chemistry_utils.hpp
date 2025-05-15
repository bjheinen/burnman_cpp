/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_CHEMISTRY_UTILS_INCLUDED
#define BURNMAN_UTILS_CHEMISTRY_UTILS_INCLUDED

#include <string>
#include <unordered_set>
#include <vector>
#include "burnman/utils/constants.hpp"

namespace utils {

  /**
   * @brief Sorts an element list to IUPAC order.
   *
   * @note Removes duplicate values.
   *
   * @param elements List of elements
   * @return ordered_elements.
   */
  inline std::vector<std::string> sort_element_list_to_IUPAC_order(
    const std::unordered_set<std::string>& unordered_elements
  ) {
    std::vector<std::string> ordered_elements;
    for (const auto& element : constants::chemistry::IUPAC_element_order) {
      if (unordered_elements.count(element)) {
        ordered_elements.push_back(element);
      }
    }
    return ordered_elements;
  }

}

#endif // BURNMAN_UTILS_CHEMISTRY_UTILS_INCLUDED
