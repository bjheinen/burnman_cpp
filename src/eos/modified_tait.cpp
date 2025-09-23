/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/eos/modified_tait.hpp"
#include <cmath>
#include <stdexcept>
#include "burnman/utils/validate_optionals.hpp"

namespace burnman::eos {

void MT::validate_parameters(types::MineralParams& params) {
  // Check for required parameters
  utils::require_set(params.V_0, "V_0");
  utils::require_set(params.K_0, "K_0");
  utils::require_set(params.Kprime_0, "Kprime_0");
  utils::require_set(params.Kdprime_0, "Kdprime_0");
  // Check for optional parameters and set defaults if missing
  utils::fallback_to_default(params.E_0, 0.0);
  utils::fallback_to_default(params.P_0, 1.0e5);
  utils::fallback_to_default(params.G_0, std::nan(""));
  utils::fallback_to_default(params.Gprime_0, std::nan(""));
  // Check values are reasonable
  utils::check_in_range(params.P_0, 0.0, 1.0e100, "P_0");
  utils::check_in_range(params.V_0, 1.0e-7, 1.0e-2, "V_0");
  utils::check_in_range(params.K_0, 1.0e9, 1.0e13, "K_0");
  utils::check_in_range(params.Kprime_0, 0.0, 40.0, "Kprime_0");
  utils::check_in_range(params.G_0, 0.0, 1.0e13, "G_0");
  utils::check_in_range(params.Gprime_0, -5.0, 10.0, "Gprime_0");
}

// Static member functions
TaitConstants MT::compute_tait_constants(const types::MineralParams& params) {
  // TODO: Should probably cache in some way
  //       but at the moment no way to track changes in params and
  //       no init for the EOS... hashmap may have more overhead
  //       than just computing the function each time.
  double K_0 = *params.K_0;
  double Kprime_0 = *params.Kprime_0;
  double Kdprime_0 = *params.Kdprime_0;
  double Kprime_0_plus_one = 1.0 + Kprime_0;
  double a = Kprime_0_plus_one / (Kprime_0_plus_one + K_0 * Kdprime_0);
  double b = Kprime_0 / K_0 - Kdprime_0 / Kprime_0_plus_one;
  double c = (Kprime_0_plus_one + K_0 * Kdprime_0)
    / (Kprime_0 * Kprime_0 + Kprime_0 - K_0 * Kdprime_0);
  return {a, b, c};
}

double MT::compute_modified_tait_pressure(
  double compression,
  const types::MineralParams& params
) {
  auto [a, b, c] = compute_tait_constants(params);
  double x_term = (compression + a - 1.0) / a;
  return (std::pow(x_term, -1.0 / c) - 1.0) / b + *params.P_0;
}

double MT::compute_modified_tait_bulk_modulus(
  double pressure,
  const types::MineralParams& params
) {
  auto [a, b, c] = compute_tait_constants(params);
  double p_term = 1.0 + b * (pressure - (*params.P_0));
  return (*params.K_0) * p_term * (a + (1.0 - a) * std::pow(p_term, c));
}

double MT::compute_modified_tait_volume(
  double pressure,
  const types::MineralParams& params
) {
  auto [a, b, c] = compute_tait_constants(params);
  double p_term = 1.0 + b * (pressure - *params.P_0);
  return (*params.V_0) * (1.0 - a * (1.0 - std::pow(p_term, -1.0*c)));
}

// Specific EOS functions
double MT::compute_volume(
  double pressure,
  double temperature [[maybe_unused]],
  const types::MineralParams& params
) const {
  return compute_modified_tait_volume(pressure, params);
}

double MT::compute_pressure(
  double temperature [[maybe_unused]],
  double volume,
  const types::MineralParams& params
) const {
  return compute_modified_tait_pressure(volume / *params.V_0, params);
}

double MT::compute_isothermal_bulk_modulus_reuss(
  double pressure,
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params
) const {
  return compute_modified_tait_bulk_modulus(pressure, params);
}

double MT::compute_isentropic_bulk_modulus_reuss(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  // TODO: maybe just return K_T??
  return 1.0e99;
}

double MT::compute_molar_internal_energy(
  double pressure,
  double temperature,
  double volume,
  const types::MineralParams& params
) const {
  return compute_gibbs_free_energy(pressure, temperature, volume, params)
    - volume * pressure;
}

double MT::compute_gibbs_free_energy(
  double pressure,
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params
) const {
  auto [a, b, c] = compute_tait_constants(params);
  double P_diff = pressure - *params.P_0;
  double intVdP = *params.V_0
    * (
      a / (b * (1.0 - c))
      * (std::pow (b * P_diff + 1.0, 1.0 - c) - 1.0)
      + (1.0 - a) * P_diff
    );
  return intVdP + (*params.E_0) + (*params.V_0) * (*params.P_0);
}

double MT::compute_shear_modulus(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double MT::compute_molar_heat_capacity_v(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 1.0e99;
}

double MT::compute_molar_heat_capacity_p(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 1.0e99;
}

double MT::compute_thermal_expansivity(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double MT::compute_entropy(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double MT::compute_grueneisen_parameter(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

} // namespace burnman::eos
