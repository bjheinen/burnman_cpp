/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_EOS_EINSTEIN_HPP_INCLUDED
#define BURNMAN_EOS_EINSTEIN_HPP_INCLUDED

#include "burnman/utils/eos.hpp"

/**
 * TODO: Add to EOS module
 * @namespace einstein
 * @brief Functions for the Einstein model of a solid.
 *
 * @note All functions assume SI units for all properties.
 */
namespace einstein {

  /**
   * Wrapper type to force explicit calls of Einstein functions with
   * double f instead of int napfu.
   */
  struct ExplicitDouble {
    double value;
    explicit ExplicitDouble(double v) : value(v) {}
  };

  /**
   * @brief Computes thermal energy of material.
   *
   * @param temperature The temperature to evaluate [K].
   * @param einstein_temperature Einstein temperature [K].
   * @param napfu Number of atoms per molecule.
   *
   * @return Thermal energy in [J/mol].
   */
  double compute_thermal_energy(
    double temperature,
    double einstein_temperature,
    int napfu);

  /**
   * @brief Computes thermal energy of material.
   *
   * @param temperature The temperature to evaluate [K].
   * @param einstein_temperature Einstein temperature [K].
   * @param f Excess correction parameter.
   *
   * @return Thermal energy in [J/mol].
   */
  double compute_thermal_energy(
    double temperature,
    double einstein_temperature,
    ExplicitDouble f);

  // Deleted overload to prevent implicit conversion
  double compute_thermal_energy(
    double temperature,
    double einstein_temperature,
    double napfu) = delete;

  // Internal implementation
  double compute_thermal_energy_impl(
    double temperature,
    double einstein_temperature,
    double napfu);

  /**
   * @brief Computes the molar heat capacity at constant volume.
   *
   * @param temperature In [K].
   * @param einstein_temperature Einstein T in [K].
   * @param napfu Number of atoms per formua unit.
   *
   * @return Heat capacity at constant volume in [J/K/mol].
   */
  double compute_molar_heat_capacity_v(
    double temperature,
    double einstein_temperature,
    int napfu);

  /**
   * @brief Computes the molar heat capacity at constant volume.
   *
   * @param temperature In [K].
   * @param einstein_temperature Einstein T in [K].
   * @param f Excess correction parameter.
   *
   * @return Heat capacity at constant volume in [J/K/mol].
   */
  double compute_molar_heat_capacity_v(
    double temperature,
    double einstein_temperature,
    ExplicitDouble f);

  // Deleted overload to prevent implicit conversion
  double compute_molar_heat_capacity_v(
    double temperature,
    double einstein_temperature,
    double napfu) = delete;

  // Internal implementation
  double compute_molar_heat_capacity_v_impl(
    double temperature,
    double einstein_temperature,
    double napfu);

  /**
   * @brief Compute the Einstein model helmholtz free energy.
   *
   * The helmholtz free energy of lattice vibrations in the Einstein model.
   * This does NOT include the zero point energy for the lattice.
   * This will cancel as long as you are calculating relative
   * differences in F.
   *
   * @param temperature [K].
   * @param einstein_temperature [K].
   * @param napfu Number of atoms in formula unit.
   *
   * @return Helmholtz energy in [J].
   */
  double compute_helmholtz_free_energy(
    double temperature,
    double einstein_temperature,
    int napfu);

  /**
   * @brief Compute the Einstein model helmholtz free energy.
   *
   * The helmholtz free energy of lattice vibrations in the Einstein model.
   * This does NOT include the zero point energy for the lattice.
   * This will cancel as long as you are calculating relative
   * differences in F.
   *
   * @param temperature [K].
   * @param einstein_temperature [K].
   * @param f Excess correction parameter.
   *
   * @return Helmholtz energy in [J].
   */
  double compute_helmholtz_free_energy(
    double temperature,
    double einstein_temperature,
    ExplicitDouble f);

  // Deleted overload to prevent implicit conversion
  double compute_helmholtz_free_energy(
    double temperature,
    double einstein_temperature,
    double napfu) = delete;

  // Internal implementation
  double compute_helmholtz_free_energy_impl(
    double temperature,
    double einstein_temperature,
    double napfu);

  /**
   * @brief Entropy due to lattice vibrations in the Einstein model.
   *
   * @param temperature [K].
   * @param einstein_temperature [K].
   * @param napfu Number of atoms in formula unit.
   *
   * @return Entropy in [J/K].
   */
  double compute_entropy(
    double temperature,
    double einstein_temperature,
    int napfu);

  /**
   * @brief Entropy due to lattice vibrations in the Einstein model.
   *
   * @param temperature [K].
   * @param einstein_temperature [K].
   * @param f Excess correction parameter.
   *
   * @return Entropy in [J/K].
   */
  double compute_entropy(
    double temperature,
    double einstein_temperature,
    ExplicitDouble f);

  // Deleted overload to prevent implicit conversion
  double compute_entropy(
    double temperature,
    double einstein_temperature,
    double napfu) = delete;

  // Internal implementation
  double compute_entropy_impl(
    double temperature,
    double einstein_temperature,
    double napfu);

  /**
   * @brief First temperature derivative of heat capacity at constant V.
   *
   * @param temperature [K].
   * @param einstein_temperature [K].
   * @param napfu Number of atoms in formula unit.
   *
   * @return dCvdT [J/K^2/mol].
   */
  double compute_dmolar_heat_capacity_v_dT(
    double temperature,
    double einstein_temperature,
    int napfu);

  /**
   * @brief First temperature derivative of heat capacity at constant V.
   *
   * @param temperature [K].
   * @param einstein_temperature [K].
   * @param f Excess correction parameter.
   *
   * @return dCvdT [J/K^2/mol].
   */
  double compute_dmolar_heat_capacity_v_dT(
    double temperature,
    double einstein_temperature,
    ExplicitDouble f);

  // Deleted overload to prevent implicit conversion
  double compute_dmolar_heat_capacity_v_dT(
    double temperature,
    double einstein_temperature,
    double napfu) = delete;

  // Internal implementation
  double compute_dmolar_heat_capacity_v_dT_impl(
    double temperature,
    double einstein_temperature,
    double napfu);

}

#endif // BURNMAN_EOS_EINSTEIN_HPP_INCLUDED
