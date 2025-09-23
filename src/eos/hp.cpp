/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/eos/hp.hpp"
#include <cmath>
#include <stdexcept>
#include "burnman/utils/validate_optionals.hpp"
#include "burnman/eos/components/einstein.hpp"
#include "burnman/eos/modified_tait.hpp"

namespace burnman::eos {

void HP_TMT::validate_parameters(types::MineralParams& params) {
  // Check for required keys
  utils::require_set(params.V_0, "V_0");
  utils::require_set(params.K_0, "K_0");
  utils::require_set(params.Kprime_0, "Kprime_0");
  utils::require_set(params.Kdprime_0, "Kdprime_0");
  utils::require_set(params.napfu, "napfu");
  utils::require_set(params.molar_mass, "molar_mass");
  utils::require_set(params.Cp, "Cp params");
  utils::require_set(params.a_0, "a_0");
  // Set defaults for missing values
  utils::fallback_to_default(params.E_0, 0.0);
  utils::fallback_to_default(params.P_0, 1.0e5);
  utils::fallback_to_default(params.T_0, 298.15);
  utils::fallback_to_default(params.H_0, std::nan(""));
  utils::fallback_to_default(params.S_0, std::nan(""));
  utils::fallback_to_default(params.G_0, std::nan(""));
  utils::fallback_to_default(params.Gprime_0, std::nan(""));
  // Approx T_e from HP2011, p.346, par.1
  if (!params.T_einstein.has_value()) {
    // TODO: watch for potential bug
    //  if using S_0 as fitting param --> will need to update T_e
    params.T_einstein = 10636.0 / (*params.S_0 / *params.napfu + 6.44);
  }
  // Check values are reasonable
  utils::check_in_range(params.T_0, 0.0, 1.0e5, "T_0");
  utils::check_in_range(params.P_0, 0.0, 1.0e100, "P_0");
  utils::check_in_range(params.V_0, 1.0e-7, 1.0e-2, "V_0");
  utils::check_in_range(params.K_0, 1.0e9, 1.0e13, "K_0");
  utils::check_in_range(params.Kprime_0, 0.0, 40.0, "Kprime_0");
  utils::check_in_range(params.G_0, 0.0, 1.0e13, "G_0");
  utils::check_in_range(params.Gprime_0, -5.0, 10.0, "Gprime_0");
  utils::check_in_range(params.S_0, 0.0, 1.0e3, "S_0");
  utils::check_in_range(params.a_0, 0.0, 1.0e-3, "a_0");
  utils::check_in_range(params.napfu, 1, 1000, "napfu");
  utils::check_in_range(params.molar_mass, 0.001, 10.0, "molar_mass");
  // TODO: checks for negative Cp_0 and Cp_2000
}

double HP_TMT::compute_volume(
  double pressure,
  double temperature,
  const types::MineralParams& params
) const {
  double Pth = compute_relative_thermal_pressure(temperature, params);
  return MT::compute_modified_tait_volume(pressure - Pth, params);
}

double HP_TMT::compute_pressure(
  double temperature,
  double volume,
  const types::MineralParams& params
) const {
  // TODO: Check requires V/V0
  double Pth = compute_relative_thermal_pressure(temperature, params);
  return MT::compute_modified_tait_pressure(volume / *params.V_0, params)
    + Pth;
}

double HP_TMT::compute_grueneisen_parameter(
  double pressure,
  double temperature,
  double volume,
  const types::MineralParams& params
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
  double volume [[maybe_unused]],
  const types::MineralParams& params
) const {
  double Pth = compute_relative_thermal_pressure(temperature, params);
  return MT::compute_modified_tait_bulk_modulus(pressure - Pth, params);
}

double HP_TMT::compute_isentropic_bulk_modulus_reuss(
  double pressure,
  double temperature,
  double volume,
  const types::MineralParams& params
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
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double HP_TMT::compute_molar_heat_capacity_v(
  double pressure,
  double temperature,
  double volume,
  const types::MineralParams& params
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
  double volume [[maybe_unused]],
  const types::MineralParams& params
) const {
  auto [a, b, c] = MT::compute_tait_constants(params);
  double Pth = compute_relative_thermal_pressure(temperature, params);
  double Cv_einst = einstein::compute_molar_heat_capacity_v(
    temperature, *params.T_einstein, *params.napfu);
  double CT_over_CTref = Cv_einst
    / einstein::compute_molar_heat_capacity_v(
      *params.T_0, *params.T_einstein, *params.napfu);
  double P_iso = pressure - *params.P_0 - Pth;
  double dintVdPdT = (*params.V_0) * (*params.a_0) * (*params.K_0) * a
    * CT_over_CTref
    * (std::pow((1.0 + b * P_iso), -c)
      - std::pow((1.0 - b * Pth), -c));
  double dSdT0 = (*params.V_0) * (*params.K_0)
    * CT_over_CTref * CT_over_CTref * (*params.a_0) * (*params.a_0)
    * (std::pow((1.0 + b * P_iso), -1.0 - c)
      - std::pow((1.0 + b * (-Pth)), -1.0 - c));
  double x = *params.T_einstein / temperature;
  double dCv_einstdT = - (Cv_einst
    * (1.0 - 2.0 / x + 2.0 / (std::exp(x) - 1.0))
    * x / temperature);
  double dSdT1 = -dintVdPdT * dCv_einstdT / Cv_einst;
  double dSdT = dSdT0 + dSdT1;
  return compute_molar_heat_capacity_pref(temperature, params)
    + temperature * dSdT;
}

double HP_TMT::compute_thermal_expansivity(
  double pressure,
  double temperature,
  double volume [[maybe_unused]],
  const types::MineralParams& params
) const {
  auto [a, b, c] = MT::compute_tait_constants(params);
  double Pth = compute_relative_thermal_pressure(temperature, params);
  double P_iso = pressure - Pth - *params.P_0;
  double C_V0 = einstein::compute_molar_heat_capacity_v(
    *params.T_0, *params.T_einstein, *params.napfu);
  double C_V = einstein::compute_molar_heat_capacity_v(
    temperature, *params.T_einstein, *params.napfu);
  double bPiso_plus_one = 1.0 + b * P_iso;
  return *params.a_0 * (C_V / C_V0)
    / (bPiso_plus_one * (a + (1.0 - a) * std::pow(bPiso_plus_one, c)));
}

double HP_TMT::compute_gibbs_free_energy(
  double pressure,
  double temperature,
  double volume [[maybe_unused]],
  const types::MineralParams& params
) const {
  auto [a, b, c] = MT::compute_tait_constants(params);
  double Pth = compute_relative_thermal_pressure(temperature, params);
  double P_iso = pressure - Pth - *params.P_0;
  double P_diff = pressure - *params.P_0;
  double intVdP;
  if (pressure != *params.P_0) {
    intVdP = P_diff * (*params.V_0)
      * (1.0 - a
        + (a
          * (std::pow((1.0 - b * Pth), 1.0 - c)
            - std::pow((1.0 + b * P_iso), 1.0 - c))
          / (b * (c - 1.0) * P_diff)));
  } else {
    intVdP = 0.0;
  }
  return *params.H_0
    + compute_intCpdT(temperature, params)
    - temperature * (*params.S_0 + compute_intCpoverTdT(temperature, params))
    + intVdP;
}

double HP_TMT::compute_entropy(
  double pressure,
  double temperature,
  double volume [[maybe_unused]],
  const types::MineralParams& params
) const {
  // S(P_0, T) = S_0 + int Cp/T dT
  // S(P, T) = S(P_0, T) + dintVdP/dT
  auto [a, b, c] = MT::compute_tait_constants(params);
  double Pth = compute_relative_thermal_pressure(temperature, params);
  double P_iso = pressure - Pth - *params.P_0;
  // C_V0(T) / C_V0(T_ref)
  double CT_over_CTref = einstein::compute_molar_heat_capacity_v(
      temperature, *params.T_einstein, *params.napfu)
    / einstein::compute_molar_heat_capacity_v(
      *params.T_0, *params.T_einstein, *params.napfu);
  double dintVdPdT = (*params.V_0) * (*params.a_0) * (*params.K_0) * a
    * CT_over_CTref
    * (std::pow((1.0 + b * P_iso), -c)
      - std::pow((1.0 - b * Pth), -c));
  return *params.S_0
    + compute_intCpoverTdT(temperature, params)
    + dintVdPdT;
}

double HP_TMT::compute_helmholtz_free_energy(
  double pressure,
  double temperature,
  double volume,
  const types::MineralParams& params
) const {
  // TODO: check if call to compute_volume needed instead
  // TODO: these are unused, dead code?
  //       --> calculated via Mineral instead.
  // Python --> uses V calculated from P instead
  // --> remove in latest version as everything done in Mineral....
  return compute_gibbs_free_energy(
      pressure, temperature, volume, params)
    - pressure * volume;
}

double HP_TMT::compute_enthalpy(
  double pressure,
  double temperature,
  double volume,
  const types::MineralParams& params
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
  const types::MineralParams& params
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
  const types::MineralParams& params
) const {
  types::CpParams cp_params = *params.Cp;
  return cp_params.a
    + cp_params.b * temperature
    + cp_params.c / (temperature * temperature)
    + cp_params.d / std::sqrt(temperature);
}

double HP_TMT::compute_intCpdT(
  double temperature,
  const types::MineralParams& params
) const {
  double T_0 = *params.T_0;
  types::CpParams cp_params = *params.Cp;
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
  const types::MineralParams& params
) const {
  types::CpParams cp_params = *params.Cp;
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
  const types::MineralParams& params
) const {
  return compute_thermal_pressure(temperature, params)
    - compute_thermal_pressure(*params.T_0, params);
}

double HP_TMT::compute_thermal_pressure(
  double temperature,
  const types::MineralParams& params
) const {
  double E_th = einstein::compute_thermal_energy(
    temperature, *params.T_einstein, *params.napfu);
  double C_V0 = einstein::compute_molar_heat_capacity_v(
    *params.T_0, *params.T_einstein, *params.napfu);
  return (*params.a_0) * (*params.K_0) / C_V0 * E_th;
}

} // namespace burnman::eos
