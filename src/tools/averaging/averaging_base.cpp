/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/tools/averaging/averaging_base.hpp"

namespace averaging {

double AveragingScheme::average_density(
  const Eigen::ArrayXd& volumes,
  const Eigen::ArrayXd& densities
) const {
  return (volumes * densities).sum() / volumes.sum();
}

double AveragingScheme::average_thermal_expansivity(
  const Eigen::ArrayXd& volumes,
  const Eigen::ArrayXd& alphas
) const {
  return (volumes * alphas).sum() / volumes.sum();
}

double AveragingScheme::average_heat_capacity_v(
  const Eigen::ArrayXd& fractions,
  const Eigen::ArrayXd& c_v
) const {
  return (fractions * c_v).sum();
}

double AveragingScheme::average_heat_capacity_p(
  const Eigen::ArrayXd& fractions,
  const Eigen::ArrayXd& c_p
) const {
  return (fractions * c_p).sum();
}

} // namespace averaging
