/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_EOS_BUKOWINSKI_ELECTRONIC_HPP_INCLUDED
#define BURNMAN_EOS_BUKOWINSKI_ELECTRONIC_HPP_INCLUDED

#include "burnman/utils/types/mineral_params.hpp"

namespace burnman {
namespace eos {

/**
 * @namespace burnman::eos::bukowinski
 * @brief Functions for the Bukowinski model for the electronic
 *        component of the Helmholtz free energy (as used by
 *        Stixrude & Lithgow-Bertelloni, 2024).
 *
 * @note All functions assume SI units for all properties.
 */
namespace bukowinski {
  /**
   * @brief Electronic component of Helmholtz free energy.
   *
   * @param temperature The temperature to evaluate [K].
   * @param volume [m^3].
   * @param params MineralParams object reference.
   *
   * @return Electronic component of Helmholtz free energy [J].
   */
  double compute_helmholtz_el(
    double temperature,
    double volume,
    const MineralParams& params);

  /**
   * @brief Electronic component of pressure, -dF/dV.
   *
   * @param temperature The temperature to evaluate [K].
   * @param volume [m^3].
   * @param params MineralParams object reference.
   *
   * @return -dF/dV [Pa].
   */
  double compute_pressure_el(
    double temperature,
    double volume,
    const MineralParams& params);

  /**
   * @brief Electronic component entropy, -dF/dT.
   *
   * @param temperature The temperature to evaluate [K].
   * @param volume [m^3].
   * @param params MineralParams object reference.
   *
   * @return -dF/dT [J/K].
   */
  double compute_entropy_el(
    double temperature,
    double volume,
    const MineralParams& params);

  /**
   * @brief K_T / V, where K_T = -V dP/dV.
   *
   * @param temperature The temperature to evaluate [K].
   * @param volume [m^3].
   * @param params MineralParams object reference.
   *
   * @return K_T/V
   */
  double compute_KT_over_V(
    double temperature,
    double volume,
    const MineralParams& params);

  /**
   * @brief C_v / T, where C_v = T dS/dT.
   *
   * @param temperature The temperature to evaluate [K].
   * @param volume [m^3].
   * @param params MineralParams object reference.
   *
   * @return dS/dT
   */
  double compute_CV_over_T(
    double volume,
    const MineralParams& params);

  /**
   * @brief alpha K_T = dP/dT.
   *
   * @param temperature The temperature to evaluate [K].
   * @param volume [m^3].
   * @param params MineralParams object reference.
   *
   * @return dP/dT
   */
  double compute_alpha_KT(
    double temperature,
    double volume,
    const MineralParams& params);

} // namespace bukowinski
} // namespace eos
} // namespace burnman

#endif // BURNMAN_EOS_BUKOWINSKI_ELECTRONIC_HPP_INCLUDED
