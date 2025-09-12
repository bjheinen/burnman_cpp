/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_TOOLS_AVERAGING_REUSS_HPP_INCLUDED
#define BURNMAN_TOOLS_AVERAGING_REUSS_HPP_INCLUDED

#include "burnman/tools/averaging/averaging_base.hpp"

namespace burnman {
namespace averaging {

/**
 * @class Reuss
 * @brief Reuss averaging scheme.
 *
 * Class for computing the Reuss (iso-stress) bound for elastic
 * properties.
 *
 * Overrides bulk and shear moduli averaging in `AveragingScheme'.
 */
class Reuss : public AveragingScheme {
 public:
  double average_bulk_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;

  double average_shear_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;
};

} // namespace averaging
} // namespace burnman

#endif // BURNMAN_TOOLS_AVERAGING_REUSS_HPP_INCLUDED
