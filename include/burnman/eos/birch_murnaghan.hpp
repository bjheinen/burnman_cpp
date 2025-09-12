/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_EOS_BIRCH_MURNAGHAN_HPP_INCLUDED
#define BURNMAN_EOS_BIRCH_MURNAGHAN_HPP_INCLUDED

#include "burnman/utils/types/mineral_params.hpp"
#include "burnman/core/equation_of_state.hpp"

namespace burnman {
namespace eos {

/**
 * @class BM3
 * @brief Base class for the Birch Murnaghan equation of state.
 *
 * The EOS is third order in strain, and has no temperature dependence.
 * Uses third order expansion also for shear modulus. To fit shear modulus
 * to second order expansion @see BM2
 *
 * @note All functions assume SI units for all properties.
 */
class BM3 : public EquationOfState{
 public:

  // Helper functions
  bool validate_parameters(types::MineralParams& params) override;

  // Static functions (for public access outside class)
  /**
   * @brief Evaluate the BM EOS pressure.
   *
   * @param inv_compression V_0/V.
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Pressure in [Pa].
   */
  static double compute_birch_murnaghan(
    double inv_compression,
    const types::MineralParams& params);

  /**
   * @brief Evaluate the bulk modulus, K
   *
   * @param volume Volume to evaluate [cm^3].
   * @param params types::MineralParams object.
   *
   * @return Bulk modulus in [Pa].
   */
  static double compute_bm_bulk_modulus(
    double volume,
    const types::MineralParams& params);

  /**
   * @brief Third order exapansion for shear modulus
   *
   * @param volume Volume to evaluate [cm^3].
   * @param params types::MineralParams object.
   *
   * @return Shear modulus in [Pa].
   */
  static double compute_third_order_shear_modulus(
    double volume,
    const types::MineralParams& params);

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

  /**
   * @copydoc EquationOfState::compute_shear_modulus
   *
   * @note Third order expansion
   */
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

  double compute_molar_internal_energy(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const override;

 private:
  /**
   * @brief GSL function wrapper to compute P(V) - P
   *
   * @param x Volume to test (passed by solver)
   * @param p Generic pointer for parameter object
   * @see `eos::gsl_params::SolverParams_P`
   */
  static double bm_gsl_wrapper(double x, void* p);

};

/**
 * @class BM2
 * @brief Birch-Murnaghan EOS with second order expansion for shear modulus.
 *
 * Derived from BM3.
 * The EOS is third order in strain, and has no temperature dependence.
 * Overrides compute_shear_modulus to use the second order expansion.
 *
 * @note All functions assume SI units for all properties.
 */
class BM2 : public BM3{
 public:
  /**
   * @brief Second order exapansion for shear modulus
   *
   * @param volume Volume to evaluate [cm^3].
   * @param params types::MineralParams object.
   *
   * @return Shear modulus in [Pa].
   */
  static double compute_second_order_shear_modulus(
    double volume,
    const types::MineralParams& params);

  /**
   * @copydoc EquationOfState::compute_shear_modulus
   *
   * @note Second order expansion
   */
  double compute_shear_modulus(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const override;
};

} // namespace eos
} // namespace burnman

#endif // BURNMAN_EOS_BIRCH_MURNAGHAN_HPP_INCLUDED