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

namespace burnman::eos {

bool MT::validate_parameters(MineralParams& params) {
  // Check for required keys
  if (!params.V_0.has_value()) {
    throw std::invalid_argument("params object missing parameter: V_0");
  }
  if (!params.K_0.has_value()) {
    throw std::invalid_argument("params object missing parameter: K_0");
  }
  if (!params.Kprime_0.has_value()) {
    throw std::invalid_argument("params object missing parameter: Kprime_0");
  }
  if (!params.Kdprime_0.has_value()) {
    throw std::invalid_argument("params object missing parameter: Kdprime_0");
  }
  // Set defaults for missing values
  if (!params.E_0.has_value()) {
    params.E_0 = 0.0;
  }
  if (!params.P_0.has_value()) {
    params.P_0 = 1.0e5;
  }
  // Set G to NaN unless user has set
  if (!params.G_0.has_value()) {
    params.G_0 = std::nan("");
  }
  if (!params.Gprime_0.has_value()) {
    params.Gprime_0 = std::nan("");
  }
  // Check values are reasonable
  // TODO: warnings?
  if (*params.P_0 < 0.0) {
    ; // warnings.warn("Unusual value for P_0", stacklevel=2)
  }
  if (*params.V_0 < 1.0e-7 || *params.V_0 > 1.0e-2) {
    ; //warning warnings.warn("Unusual value for V_0", stacklevel=2)
  }
  if (*params.K_0 < 1.0e9 || *params.K_0 > 1.0e13) {
    ; // warning warnings.warn("Unusual value for K_0", stacklevel=2)
  }
  if (*params.Kprime_0 < 0.0 || *params.Kprime_0 > 40.0) {
    ; // warning warnings.warn("Unusual value for Kprime_0", stacklevel=2)
  }
  if (*params.G_0 < 0.0 || *params.G_0 > 1.0e13) {
    ; // warnings.warn("Unusual value for G_0", stacklevel=2)
  }
  if (*params.Gprime_0 < -5.0 || *params.Gprime_0 > 10.0) {
    ; // warnings.warn("Unusual value for Gprime_0", stacklevel=2)
  }
  return 1;
}

// Static member functions
TaitConstants MT::compute_tait_constants(const MineralParams& params) {
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
  const MineralParams& params
) {
  auto [a, b, c] = compute_tait_constants(params);
  double x_term = (compression + a - 1.0) / a;
  return (std::pow(x_term, -1.0 / c) - 1.0) / b + *params.P_0;
}

double MT::compute_modified_tait_bulk_modulus(
  double pressure,
  const MineralParams& params
) {
  auto [a, b, c] = compute_tait_constants(params);
  double p_term = 1.0 + b * (pressure - (*params.P_0));
  return (*params.K_0) * p_term * (a + (1.0 - a) * std::pow(p_term, c));
}

double MT::compute_modified_tait_volume(
  double pressure,
  const MineralParams& params
) {
  auto [a, b, c] = compute_tait_constants(params);
  double p_term = 1.0 + b * (pressure - *params.P_0);
  return (*params.V_0) * (1.0 - a * (1.0 - std::pow(p_term, -1.0*c)));
}

// Specific EOS functions
double MT::compute_volume(
  double pressure,
  double temperature [[maybe_unused]],
  const MineralParams& params
) const {
  return compute_modified_tait_volume(pressure, params);
}

double MT::compute_pressure(
  double temperature [[maybe_unused]],
  double volume,
  const MineralParams& params
) const {
  return compute_modified_tait_pressure(volume / *params.V_0, params);
}

double MT::compute_isothermal_bulk_modulus_reuss(
  double pressure,
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params
) const {
  return compute_modified_tait_bulk_modulus(pressure, params);
}

double MT::compute_isentropic_bulk_modulus_reuss(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params [[maybe_unused]]
) const {
  // TODO: maybe just return K_T??
  return 1.0e99;
}

double MT::compute_molar_internal_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return compute_gibbs_free_energy(pressure, temperature, volume, params)
    - volume * pressure;
}

double MT::compute_gibbs_free_energy(
  double pressure,
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params
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
  const MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double MT::compute_molar_heat_capacity_v(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params [[maybe_unused]]
) const {
  return 1.0e99;
}

double MT::compute_molar_heat_capacity_p(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params [[maybe_unused]]
) const {
  return 1.0e99;
}

double MT::compute_thermal_expansivity(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double MT::compute_entropy(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double MT::compute_grueneisen_parameter(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

} // namespace burnman::eos
