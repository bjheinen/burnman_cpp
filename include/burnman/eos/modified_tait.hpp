/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_EOS_MODIFIED_TAIT_HPP_INCLUDED
#define BURNMAN_EOS_MODIFIED_TAIT_HPP_INCLUDED

#include "burnman/utils/types/mineral_params.hpp"
#include "burnman/core/equation_of_state.hpp"

namespace burnman {
namespace eos {

/**
 * @brief Struct to hold tait constants
 */
struct TaitConstants {
  double a, b, c;
};

/**
* @class MT
* @brief Base class for the generic modified tait equation of state.
*
* References for the EOS can be found in Huang & Chow, 1974 and
* Holland & Powell, 2011.
*
* @note All functions assume SI units for all properties.
*/
class MT : public EquationOfState{
 public:

  // Helper functions
  bool validate_parameters(MineralParams& params) override;

  // Static functions (for public access outside class)
  /**
  * @brief Evaluate the MT EOS pressure.
  *
  * @param compression V/V_0.
  * @param params Mineral parameters object of type MineralParams
  *
  * @return Pressure in [Pa].
  */
  static double compute_modified_tait_pressure(
    double compression,
    const MineralParams& params);

  /**
  * @brief Evaluate the isothermal bulk modulus, K_T.
  *
  * @param pressure Pressure to evaluate [Pa].
  * @param params MineralParams object.
  *
  * @return Bulk modulus in [Pa].
  */
  static double compute_modified_tait_bulk_modulus(
    double pressure,
    const MineralParams& params);

  /**
  * @brief Calculates volume from pressure.
  *
  * @param pressure Pressure to evaluate [Pa].
  * @param params MineralParams object.
  *
  * @return Volume in [m^3].
  */
  static double compute_modified_tait_volume(
    double volume,
    const MineralParams& params);

  /**
   * @brief Expands MT parameters a, b, c.
   *
   * a, b, c are derived from K_T and its two first pressure
   * derivatives. Eq. 4 in Holland & Powell, 2011.
   *
   * @param params MineralParams object ref.
   *
   * @return Struct holding a, b, c
   */
  static TaitConstants compute_tait_constants(const MineralParams& params);

  // Specific EOS functions
  double compute_volume(
    double pressure,
    double temperature,
    const MineralParams& params) const override;

  double compute_pressure(
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_isothermal_bulk_modulus_reuss(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_isentropic_bulk_modulus_reuss(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_molar_internal_energy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_gibbs_free_energy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_shear_modulus(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_molar_heat_capacity_v(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_molar_heat_capacity_p(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_thermal_expansivity(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_entropy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_grueneisen_parameter(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

};

} // namespace eos
} // namespace burnman

#endif // BURNMAN_EOS_MODIFIED_TAIT_HPP_INCLUDED
