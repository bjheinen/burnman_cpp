/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/tools/averaging/reuss.hpp"
#include "burnman/tools/averaging/averaging_utils.hpp"

namespace averaging {

double Reuss::average_bulk_moduli(
  const Eigen::ArrayXd& volumes,
  const Eigen::ArrayXd& bulk_moduli,
  const Eigen::ArrayXd& shear_moduli [[maybe_unused]]
) const {
  return utils::reuss_fn(volumes, bulk_moduli);
}

double Reuss::average_shear_moduli(
  const Eigen::ArrayXd& volumes,
  const Eigen::ArrayXd& bulk_moduli [[maybe_unused]],
  const Eigen::ArrayXd& shear_moduli
) const {
  return utils::reuss_fn(volumes, shear_moduli);
}

} // namespace averaging
