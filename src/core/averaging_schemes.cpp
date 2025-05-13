/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */

#include "burnman/core/averaging_schemes.hpp"
#include "burnman/utils/constants.hpp"

double averaging::reuss_fn(
  const Eigen::ArrayXd& phase_volumes,
  const Eigen::ArrayXd& X
) {
  Eigen::ArrayXd vol_frac = phase_volumes / phase_volumes.sum();
  double eps = constants::precision::abs_tolerance;
  // Warn when X <= 0 and |vol_frac| >= eps 
  Eigen::Array<bool, Eigen::Dynamic, 1> problematic =
      (X <= 0.0).select(vol_frac.abs() > eps, false);
  if (problematic.any()) {
      // Warn!
      // Could use cerr?
      // std::cerr << "Oops, called reuss_average with Xi <= 0!" << std::endl;
      return 0.0;
  } else {
    return 1.0 / (vol_frac / X).sum();
  }
}
