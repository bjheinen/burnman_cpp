/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_AVERAGING_SCHEMES_HPP_INCLUDED
#define BURNMAN_CORE_AVERAGING_SCHEMES_HPP_INCLUDED

#include <Eigen/Dense>

/**
 * TODO: Add to module?
 * @namespace averaging
 * @brief Functions for computing averages.
 *
 * TODO: Could move to utils?
 * Function for use in averaging schemes, VRH, etc.
 *
 */
namespace averaging {

  /**
   * @brief Calculated Reuss (iso-stress) average.
   *
   * @param phase_volumes Partial volumes or molar fractions
   * @param X Phase moduli to average
   *
   * @return Reuss average
   */

  double reuss_fn(
    const Eigen::ArrayXd& phase_volumes,
    const Eigen::ArrayXd& X);

} // namespace averaging

#endif // BURNMAN_CORE_AVERAGING_SCHEMES_HPP_INCLUDED
