/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_EOS_MIE_GRUENEISEN_DEBYE_HPP_INCLUDED
#define BURNMAN_EOS_MIE_GRUENEISEN_DEBYE_HPP_INCLUDED

#include "burnman/utils/eos.hpp"
#include "burnman/core/equation_of_state.hpp"

/**
 * @class MGD3
 * @brief Class for the Mie-Grueneisen-Debye equation of state.
 *
 * Uses finite-strain Birch-Murnaghan for the isothermal portion.
 * References can be found in many places, e.g.
 *   Shim, Duffy and Kenichi (2002)
 *   Jackson and Rigden (1996)
 *   Matas et al. (2007) Appendices
 * The thermal correction to the shear modulus is from Hama and Suito (1998).
 * Uses third order expansion also for shear modulus. To fit shear modulus
 * to second order expansion @see MGD2
 * 
 * @note All functions assume SI units for all properties.
 */
class MGD3 : public EquationOfState{
 public:

  // Helper functions
  bool validate_parameters(MineralParams& params) override;

  // Specific EOS functions
  /**
   * @copydoc EquationOfState::compute_volume
   *
   * @note Matas et al. eq. B7.
   */
  double compute_volume(
    double pressure,
    double temperature,
    const MineralParams& params) const override;

  /**
   * @copydoc EquationOfState::compute_pressure
   *
   * @note Matas et al. eq. B7.
   */
  double compute_pressure(
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_grueneisen_parameter(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  /**
   * @copydoc EquationOfState::compute_isothermal_bulk_modulus_reuss
   *
   * @note Matas et al. eq. B8, B13.
   */
  double compute_isothermal_bulk_modulus_reuss(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  /**
   * @copydoc EquationOfState::compute_isentropic_bulk_modulus_reuss
   *
   * @note Matas et al. eq. D6.
   */
  double compute_isentropic_bulk_modulus_reuss(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  /**
   * @copydoc EquationOfState::compute_shear_modulus
   *
   * @note Third order expansion, Matas et al. eq. B11.
   */
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

  double compute_gibbs_free_energy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_entropy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_molar_internal_energy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_helmholtz_free_energy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_enthalpy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

 protected:

  /**
   * @brief Thermal correction to shear modulus
   *
   * @param temperature Temperature to evaluate [K].
   * @param volume Volume to evaluate [m^3].
   * @param params Mineral parameters object of type MineralParams.
   *
   * @return G_thermal [Pa].
   */
  double compute_thermal_shear_modulus(
    double temperature,
    double volume,
    const MineralParams& params) const;

 private:
  /**
   * @brief Compute the thermal correction to the bulk modulus.
   *
   * @param temperature Temperature to evaluate [K].
   * @param volume Volume to evaluate [m^3].
   * @param params Mineral parameters object of type MineralParams.
   *
   * @return K_th in [Pa].
   */
  double compute_thermal_bulk_modulus(
    double temperature,
    double volume,
    const MineralParams& params) const;

  /**
   * @brief Computes the Grueneisen parameter.
   *
   * Matas et al. eq. B6.
   *
   * @param x V_0/V.
   * @param params Mineral parameters object of type MineralParams
   * 
   * @return Grueneisen parameter [unitless].
   */
  static double compute_mgd_grueneisen_parameter(
    double x,
    const MineralParams& params);

  /**
   * @brief Compute the Debye temperature.
   *
   * See Matas et al. eq. B6.
   *
   * @param x V_0/V (inverse compression).
   * @param params Mineral parameters object of type MineralParams.
   *
   * @return Debye temperature in [K].
   */
  static double compute_debye_temperature(
    double x,
    const MineralParams& params);

  /**
   * @brief Compute the isotropic thermal pressure.
   *
   * See Matas et al. eq. B4.
   *
   * @param temperature Temperature to evaluate [K].
   * @param volume Volume to evaluate [m^3].
   * @param params Mineral parameters object of type MineralParams.
   *
   * @return Thermal pressure in [Pa].
   */
  static double compute_thermal_pressure(
    double temperature,
    double volume,
    const MineralParams& params);

  /**
   * @brief GSL function wrapper to compute P(V) - P
   * 
   * @param x Volume to test (passed by solver)
   * @param p Generic pointer for parameter object
   * @see `ParamsGSL::SolverParams_P`
   */
  static double mgd_gsl_wrapper(double x, void* p);

};


/**
 * @class MGD2
 * @brief MGD EOS with second order expansion of shear modulus.
 *
 * Derived from MGD3.
 * Overrides compute_shear_modulus to use the second order expansion.
 *
 * @note All functions assume SI units for all properties.
 */
class MGD2 : public MGD3{
 public:
  /**
   * @copydoc EquationOfState::compute_shear_modulus
   *
   * @note Second order expansion
   */
  double compute_shear_modulus(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;
};

#endif // BURNMAN_EOS_MIE_GRUENEISEN_DEBYE_HPP_INCLUDED