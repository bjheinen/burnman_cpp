/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_TOOLS_AVERAGING_AVERAGING_BASE_HPP_INCLUDED
#define BURNMAN_TOOLS_AVERAGING_AVERAGING_BASE_HPP_INCLUDED

#include <Eigen/Dense>

namespace averaging {

/**
 * @class AveragingScheme
 * @brief Base interface for averaging schemes.
 *
 * Functions return scalar values for a list of properties and
 * volume fractions.

 * Derived classes should implement averaging of elastic moduli.
 */
class AveragingScheme {

 public:

  virtual ~AveragingScheme() = default;

  // Moduli averages need to be implemented in subclasses.
  /**
   * @brief Average bulk modulus, `K', for a composite.
   * @param volumes Volumes of each phase in [m^3].
   * @param bulk_moduli Bulk moduli of each phase in [Pa].
   * @param shear_moduli Shear moduli of each phase in [Pa].
   * @return Average bulk modulus, K, in [Pa].
   */
  virtual double average_bulk_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const = 0;

  /**
   * @brief Average shear modulus, `G', for a composite.
   * @param volumes Volumes of each phase in [m^3].
   * @param bulk_moduli Bulk moduli of each phase in [Pa].
   * @param shear_moduli Shear moduli of each phase in [Pa].
   * @return Average shear modulus, G, in [Pa].
   */
  virtual double average_shear_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const = 0;

  // Common default implementations
  /**
   * @brief Density average of a composite.
   *
   * Average assuming: \f$ \rho = \frac{\sum_i \rho_i V_i}{\sum_i V_i} \f$.
   *
   * @param volumes Volumes of each phase in [m^3].
   * @param densities Densities of each phase in [kg/m^3].
   * @return \f$ \rho \f$ in [kg/m^3].
    */
  virtual double average_density(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& densities) const;

  /**
   * @brief Average thermal expansion coefficient of a composite.
   * @param volumes Volumes of each phase in [m^3].
   * @param alphas Thermal expansivity of each phase in [1/K].
   * @return \f$ \alpha \f$ in [1/K].
    */
  virtual double average_thermal_expansivity(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& alphas) const;

  // TODO: currently redunant as C_v and C_p identical
  /**
   * @brief Average C_v of a composite.
   *
   * Averages C_v by molar fractions as in Ita, et al., 1992 [eq. 16].

   * @param fractions Molar fractions of each phase.
   * @param c_v Heat capacities at constant volume of each phase in [J/K/mol].
   * @return Average C_v in [J/K/mol].
    */
  virtual double average_heat_capacity_v(
    const Eigen::ArrayXd& fractions,
    const Eigen::ArrayXd& c_v) const;

  /**
   * @brief Average C_p of a composite.
   * @param fractions Molar fractions of each phase.
   * @param c_p Heat capacities at constant volume of each phase in [J/K/mol].
   * @return Average C_p in [J/K/mol].
    */
  virtual double average_heat_capacity_p(
    const Eigen::ArrayXd& fractions,
    const Eigen::ArrayXd& c_p) const;

};

} // namespace averaging

#endif // BURNMAN_TOOLS_AVERAGING_AVERAGING_BASE_HPP_INCLUDED
