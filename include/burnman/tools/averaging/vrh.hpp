/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_TOOLS_AVERAGING_VRH_HPP_INCLUDED
#define BURNMAN_TOOLS_AVERAGING_VRH_HPP_INCLUDED

#include "burnman/tools/averaging/averaging_base.hpp"

namespace burnman {
namespace averaging {

/**
 * @class VoigtReussHill
 * @brief VRH averaging scheme.
 *
 * Class for computing the Voigt-Reuss-Hill average for elastic properties.
 * Defined as arithmetic mean of Voigt and Reuss bounds.
 * Overrides bulk and shear moduli averaging in `AveragingScheme'.
 */
class VoigtReussHill : public AveragingScheme {
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

#endif // BURNMAN_TOOLS_AVERAGING_VRH_HPP_INCLUDED
