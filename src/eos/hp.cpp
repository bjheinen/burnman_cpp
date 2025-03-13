/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <cmath>
#include "burnman/eos/hp.hpp"
#include "burnman/eos/einstein.hpp"

double HP_TMT::compute_intCpdT(
  double temperature,
  const MineralParams& params
) const {
  double T_0 = *params.T_0;
  CpParams cp_params = *params.Cp;
  return
    (
      cp_params.a * temperature
      + 0.5 * cp_params.b * temperature * temperature
      - cp_params.c / temperature
      + 2.0 * cp_params.d * std::sqrt(temperature))
    - (
      cp_params.a * T_0
      + 0.5 * cp_params.b * T_0 * T_0
      - cp_params.c / T_0
      + 2.0 * cp_params.d * std::sqrt(T_0));
}

double HP_TMT::compute_intCpoverTdT(
  double temperature,
  const MineralParams& params
) const {
  CpParams cp_params = *params.Cp;
  double T_0 = *params.T_0;
  return 
    (
      cp_params.a * std::log(temperature)
      + cp_params.b * temperature
      - 0.5 * cp_params.c / (temperature * temperature)
      - 2.0 * cp_params.d / std::sqrt(temperature))
    - (
      cp_params.a * std::log(T_0)
      + cp_params.b * T_0
      - 0.5 * cp_params.c / (T_0 * T_0)
      - 2.0 * cp_params.d / std::sqrt(T_0));
}

double HP_TMT::compute_relative_thermal_pressure(
  double temperature,
  const MineralParams& params
) const {
  return compute_thermal_pressure(temperature, params)
    - compute_thermal_pressure(*params.T_0, params);
}

double HP_TMT::compute_thermal_pressure(
  double temperature,
  const MineralParams& params
) const {
  double E_th = einstein::compute_thermal_energy(
    temperature, *params.T_einstein, *params.napfu);
  double C_V0 = einstein::compute_molar_heat_capacity_v(
    *params.T_0, *params.T_einstein, *params.napfu);
  return (*params.a_0) * (*params.K_0) / C_V0 * E_th;
}
