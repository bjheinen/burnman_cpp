/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_EOS_HP_HPP_INCLUDED
#define BURNMAN_EOS_HP_HPP_INCLUDED

#include "burnman/utils/types/mineral_params.hpp"
#include "burnman/core/equation_of_state.hpp"

namespace burnman {
namespace eos {

// TODO: Maybe HP11 and HP11Mod instead of HP_TMT?
/**
 * @class HP_TMT
 * @brief Class for Holland & Powell's thermal EOS correction to the generic
 *        modifed Tait EOS (see `MT`). Desribed in Holland & Powell, 2011.
 *
 * TODO: Note here on P0 & calculation of Cp --> problems at high T
 *       Note to @see modified HP class.
 *
 * @note All functions assume SI units for all properties.
 */
class HP_TMT : public EquationOfState{
 public:

  // Helper functions
  void validate_parameters(types::MineralParams& params) override;

  // Specific EOS functions
  double compute_volume(
    double pressure,
    double temperature,
    const types::MineralParams& params) const override;

  double compute_pressure(
    double temperature,
    double volume,
    const types::MineralParams& params) const override;

  double compute_grueneisen_parameter(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const override;

  double compute_isothermal_bulk_modulus_reuss(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const override;

  double compute_isentropic_bulk_modulus_reuss(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const override;

  double compute_shear_modulus(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const override;

  double compute_molar_heat_capacity_v(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const override;

  double compute_molar_heat_capacity_p(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const override;

  double compute_thermal_expansivity(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const override;

  double compute_gibbs_free_energy(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const override;

  double compute_entropy(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const override;

  double compute_helmholtz_free_energy(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const override;

  double compute_enthalpy(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const override;

  // Additional functions
  /**
   * @brief C_p from C_v and Einstein model.
   *
   * Molar heat capacity at constant volume calculated using an
   * Einstein model.
   *
   * @note Only for comparison with internally self-consisten C_p.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Heat capacity at constant pressure in [J/K/mol].
   */
  double compute_molar_heat_capacity_p_einstein(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const;

 protected:

  /**
   * @brief Heat capacity at reference pressure as function of T.
   *
   * Cp = a + bT + cT^-2 + dT^-0.5 (Holland & Powell, 2011).
   *
   * @note In `HP_TMT' P_ref = P_0, for derived classes this may not be
   *       the case. e.g. @see TODO `NAME'
   *
   * @param temperature The temperature to evaluate [K].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Heat capacity at reference pressure in [J/K/mol].
   */
  double compute_molar_heat_capacity_pref(
    double temperature,
    const types::MineralParams& params) const;

  /**
   * @brief Thermal addition to standard state enthalpy at ambient P.
   *
   * @param temperature Temperature to evaluate [K].
   * @param params Mineral parameters object.
   *
   * @return [J/mol].
   */
  virtual double compute_intCpdT(
    double temperature,
    const types::MineralParams& params) const;

  /**
   * @brief Thermal addition to standard state entropy at ambient P.
   *
   * @param temperature Temperature to evaluate [K].
   * @param params Mineral parameters object.
   *
   * @return [J/K/mol].
   */
  virtual double compute_intCpoverTdT(
    double temperature,
    const types::MineralParams& params) const;

  /**
   * @brief Relative thermal pressure as function of T - T_0.
   *
   * Eq. 12 - 1 in Holland & Powell, 2011.
   *
   * @param temperature Temperature to evaluate [K].
   * @param params Mineral parameters object.
   *
   * @return Relative thermal pressure [Pa].
   */
  double compute_relative_thermal_pressure(
    double temperature,
    const types::MineralParams& params) const;

 private:

  /**
   * @brief Computes thermal pressure as a function of T.
   *
   * @param temperature [K].
   * @param params types::MineralParams object.
   *
   * @return Thermal pressure [Pa].
   */
  double compute_thermal_pressure(
    double temperature,
    const types::MineralParams& params) const;

};

} // namespace eos
} // namespace burnman

#endif // BURNMAN_EOS_HP_HPP_INCLUDED
