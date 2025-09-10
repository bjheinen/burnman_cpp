/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_EOS_DEBYE_HPP_INCLUDED
#define BURNMAN_EOS_DEBYE_HPP_INCLUDED

#include "burnman/utils/types/simple_types.hpp"

/**
 * TODO: Add to EOS module
 * @namespace debye
 * @brief Functions for Debye model.
 *
 * Functions required for Debye model. Used with Mie-Grueneisen and
 * Birch-Murnaghan for a full EOS.
 *
 * @note All functions assume SI units for all properties.
 */
namespace debye {

  /**
   @brief Integrand of third-order Debye function.
   */
  double debye_fn_integrand(double xi, void*); // TODO - internal use only - maybe put in unnamed namespace in .cpp only

  /**
   * @brief Evaluates the Debye function using numerical integration.
   */
  double debye_fn_quad(double x);

  /**
   * @brief Evaluates the Debye function using a Chebyshev series
   *        expansion coupled with asymptotic solutions of the function.
   *
   * @note Uses the GSL implementation adapted in PyBurnman as debye_fn_cheb
   */
  double debye_fn_cheb(double x);

  /**
   * @brief Computes thermal energy of material.
   *
   * @param temperature The temperature to evaluate [K].
   * @param debye_temperature Debye temperature [K].
   * @param napfu Number of atoms per molecule.
   *
   * @return Thermal energy in [J/mol].
   */
  double compute_thermal_energy(
    double temperature,
    double debye_temperature,
    int napfu);

  /**
   * @brief Computes thermal energy of material.
   *
   * @param temperature The temperature to evaluate [K].
   * @param debye_temperature Debye temperature [K].
   * @param f Excess correction parameter.
   *
   * @return Thermal energy in [J/mol].
   */
  double compute_thermal_energy(
    double temperature,
    double debye_temperature,
    ExplicitDouble f);

  // Deleted overload to prevent implicit conversion
  double compute_thermal_energy(
    double temperature,
    double debye_temperature,
    double napfu) = delete;

  // Internal implementation
  double compute_thermal_energy_impl(
    double temperature,
    double debye_temperature,
    double napfu);

  /**
   * @brief Computes the molar heat capacity at constant volume.
   *
   * @param temperature In [K].
   * @param debye_temperature Debye T in [K].
   * @param napfu Number of atoms per formua unit.
   *
   * @return Heat capacity at constant volume in [J/K/mol].
   */
  double compute_molar_heat_capacity_v(
    double temperature,
    double debye_temperature,
    int napfu);

  /**
   * @brief Computes the molar heat capacity at constant volume.
   *
   * @param temperature In [K].
   * @param debye_temperature Debye T in [K].
   * @param f Excess correction parameter.
   *
   * @return Heat capacity at constant volume in [J/K/mol].
   */
  double compute_molar_heat_capacity_v(
    double temperature,
    double debye_temperature,
    ExplicitDouble f);

  // Deleted overload to prevent implicit conversion
  double compute_molar_heat_capacity_v(
    double temperature,
    double debye_temperature,
    double napfu) = delete;

  // Internal implementation
  double compute_molar_heat_capacity_v_impl(
    double temperature,
    double debye_temperature,
    double napfu);

  /**
   * @brief Compute the Debye model helmholtz free energy.
   *
   * The helmholtz free energy of lattice vibrations in the Debye model.
   * This does NOT include the zero point energy for the lattice.
   * This will cancel as long as you are calculating relative
   * differences in F.
   *
   * @param temperature [K].
   * @param debye_temperature [K].
   * @param napfu Number of atoms in formula unit.
   *
   * @return Helmholtz energy in [J].
   */
  double compute_helmholtz_free_energy(
    double temperature,
    double debye_temperature,
    int napfu);

  /**
   * @brief Compute the Debye model helmholtz free energy.
   *
   * The helmholtz free energy of lattice vibrations in the Debye model.
   * This does NOT include the zero point energy for the lattice.
   * This will cancel as long as you are calculating relative
   * differences in F.
   *
   * @param temperature [K].
   * @param debye_temperature [K].
   * @param f Excess correction parameter.
   *
   * @return Helmholtz energy in [J].
   */
  double compute_helmholtz_free_energy(
    double temperature,
    double debye_temperature,
    ExplicitDouble f);

  // Deleted overload to prevent implicit conversion
  double compute_helmholtz_free_energy(
    double temperature,
    double debye_temperature,
    double napfu) = delete;

  // Internal implementation
  double compute_helmholtz_free_energy_impl(
    double temperature,
    double debye_temperature,
    double napfu);

  /**
   * @brief Entropy due to lattice vibrations in the Debye model.
   *
   * @param temperature [K].
   * @param debye_temperature [K].
   * @param napfu Number of atoms in formula unit.
   *
   * @return Entropy in [J/K].
   */
  double compute_entropy(
    double temperature,
    double debye_temperature,
    int napfu);

  /**
   * @brief Entropy due to lattice vibrations in the Debye model.
   *
   * @param temperature [K].
   * @param debye_temperature [K].
   * @param f Excess correction parameter.
   *
   * @return Entropy in [J/K].
   */
  double compute_entropy(
    double temperature,
    double debye_temperature,
    ExplicitDouble f);

  // Deleted overload to prevent implicit conversion
  double compute_entropy(
    double temperature,
    double debye_temperature,
    double napfu) = delete;

  // Internal implementation
  double compute_entropy_impl(
    double temperature,
    double debye_temperature,
    double napfu);

  /**
   * @brief First temperature derivative of heat capacity at constant V.
   *
   * @param temperature [K].
   * @param debye_temperature [K].
   * @param napfu Number of atoms in formula unit.
   *
   * @return dCvdT [J/K^2/mol].
   */
  double compute_dmolar_heat_capacity_v_dT(
    double temperature,
    double debye_temperature,
    int napfu);

  /**
   * @brief First temperature derivative of heat capacity at constant V.
   *
   * @param temperature [K].
   * @param debye_temperature [K].
   * @param f Excess correction parameter.
   *
   * @return dCvdT [J/K^2/mol].
   */
  double compute_dmolar_heat_capacity_v_dT(
    double temperature,
    double debye_temperature,
    ExplicitDouble f);

  // Deleted overload to prevent implicit conversion
  double compute_dmolar_heat_capacity_v_dT(
    double temperature,
    double debye_temperature,
    double napfu) = delete;

  // Internal implementation
  double compute_dmolar_heat_capacity_v_dT_impl(
    double temperature,
    double debye_temperature,
    double napfu);

}

#endif // BURNMAN_EOS_DEBYE_HPP_INCLUDED
