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

#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include <Eigen/Dense>
#include "burnman/utils/constants.hpp"
#include "burnman/utils/types/simple_types.hpp"

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

  /** @brief Calculates (weighted) sum of chemical formulae.
   *
   * @param formulae Vector of chemical formulae with elements as FormulaMap.
   * @param weights Optional weights - vector of equal length to formulae.
   *
   * @return Summed formula as FormulaMap.
   */
  inline FormulaMap sum_formulae(
    const std::vector<FormulaMap>& formulae,
    const Eigen::ArrayXd& weights
  ) {
    std::size_t n = formulae.size();
    if (static_cast<std::size_t>(weights.size()) != n) {
      throw std::invalid_argument(
        "Weights length must be equal to number of formulae");
    }
    FormulaMap summed_formula;
    for (std::size_t i = 0; i < n; ++i) {
      summed_formula += formulae[i] * weights[i];
    }
    return summed_formula;
  }

  /**
   * @copydoc sum_formulae(
   *  const std::vector<FormulaMap>& formulae,
   *  const Eigen::ArrayXd& weights)
   * @overload
   */
  inline FormulaMap sum_formulae(
    const std::vector<FormulaMap>& formulae
  ) {
    Eigen::ArrayXd ones = Eigen::ArrayXd::Ones(static_cast<Eigen::Index>(formulae.size()));
    return sum_formulae(formulae, ones);
  }

}

#endif // BURNMAN_UTILS_CHEMISTRY_UTILS_INCLUDED
