/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/tools/averaging/hs_lower.hpp"
#include "burnman/tools/averaging/averaging_utils.hpp"

namespace burnman::averaging {

double HashinShtrikmanLower::average_bulk_moduli(
  const Eigen::ArrayXd& volumes,
  const Eigen::ArrayXd& bulk_moduli,
  const Eigen::ArrayXd& shear_moduli
) const {
  return utils::lower_hs_bulk_fn(volumes, bulk_moduli, shear_moduli);
}

double HashinShtrikmanLower::average_shear_moduli(
  const Eigen::ArrayXd& volumes,
  const Eigen::ArrayXd& bulk_moduli,
  const Eigen::ArrayXd& shear_moduli
) const {
  return utils::lower_hs_shear_fn(volumes, bulk_moduli, shear_moduli);
}

} // namespace burnman::averaging
