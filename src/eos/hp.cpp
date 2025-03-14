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
#include "burnman/eos/modified_tait.hpp"
#include "burnman/eos/einstein.hpp"

bool HP_TMT::validate_parameters(MineralParams& params) {

}

double HP_TMT::compute_volume(
  double pressure,
  double temperature,
  const MineralParams& params
) const {
  double Pth = compute_relative_thermal_pressure(temperature, params);
  return MT::compute_modified_tait_volume(pressure - Pth, params);
}

double HP_TMT::compute_pressure(
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double Pth = compute_relative_thermal_pressure(temperature, params);
  return MT::compute_modified_tait_pressure(volume / *params.V_0, params)
    + Pth;
}

double HP_TMT::compute_grueneisen_parameter(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double alpha = compute_thermal_expansivity(
    pressure, temperature, volume, params);
  double K_T = compute_isothermal_bulk_modulus_reuss(
    pressure, temperature, volume, params);
  double C_v = compute_molar_heat_capacity_v(
    pressure, temperature, volume, params);
  return alpha * K_T * volume / C_v;
}

double HP_TMT::compute_isothermal_bulk_modulus_reuss(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double Pth = compute_relative_thermal_pressure(temperature, params);
  return MT::compute_modified_tait_bulk_modulus(pressure - Pth, params);
}

double HP_TMT::compute_isentropic_bulk_modulus_reuss(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double K_T = compute_isothermal_bulk_modulus_reuss(
    pressure, temperature, volume, params);
  double C_p = compute_molar_heat_capacity_p(
    pressure, temperature, volume, params);
  double C_v = compute_molar_heat_capacity_v(
    pressure, temperature, volume, params);
  return K_T * C_p / C_v;
}

double HP_TMT::compute_shear_modulus(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return 0.0;
}

double HP_TMT::compute_molar_heat_capacity_v(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double C_p = compute_molar_heat_capacity_p(
    pressure, temperature, volume, params);
  double alpha = compute_thermal_expansivity(
    pressure, temperature, volume, params);
  double K_T = compute_isothermal_bulk_modulus_reuss(
    pressure, temperature, volume, params);
  return C_p - volume * temperature * alpha * alpha * K_T;
}

double HP_TMT::compute_molar_heat_capacity_p(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {

}

double HP_TMT::compute_thermal_expansivity(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  auto [a, b, c] = MT::compute_tait_constants(params);
  double Pth = compute_relative_thermal_pressure(temperature, params);
  double P_diff = pressure - Pth - *params.P_0;
  double C_V0 = einstein::compute_molar_heat_capacity_v(
    *params.T_0, *params.T_einstein, *params.napfu);
  double C_V = einstein::compute_molar_heat_capacity_v(
    temperature, *params.T_einstein, *params.napfu);
  double bPdiff_plus_one = 1.0 + b * P_diff;
  return *params.a_0 * (C_V / C_V0)
    / (bPdiff_plus_one * (a + (1.0 - a) * std::pow(bPdiff_plus_one, c)));
}

double HP_TMT::compute_gibbs_free_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {

}

double HP_TMT::compute_entropy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {

}

double HP_TMT::compute_helmholtz_free_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {

}

double HP_TMT::compute_enthalpy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double G = compute_gibbs_free_energy(
    pressure, temperature, volume, params);
  double S = compute_entropy(
    pressure, temperature, volume, params);
  return G + temperature * S;
}

double HP_TMT::compute_molar_heat_capacity_p_einstein(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double alpha = compute_thermal_expansivity(
    pressure, temperature, volume, params);
  double gamma = compute_grueneisen_parameter(
    pressure, temperature, volume, params);
  double C_v = compute_molar_heat_capacity_v(
    pressure, temperature, volume, params);
  return C_v * (1.0 + gamma * alpha * temperature);
}

double HP_TMT::compute_molar_heat_capacity_pref(
  double temperature,
  const MineralParams& params
) const {
  CpParams cp_params = *params.Cp;
  return cp_params.a
    + cp_params.b * temperature
    + cp_params.c / (temperature * temperature)
    + cp_params.d / std::sqrt(temperature);
}

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
