/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/eos/vinet.hpp"
#include <cmath>
#include <stdexcept>
#include "burnman/eos/components/gsl_params.hpp"
#include "burnman/utils/validate_optionals.hpp"
#include "burnman/optim/roots/brent.hpp"

namespace burnman::eos {

void Vinet::validate_parameters(types::MineralParams& params) {
  // Check for required parameters
  utils::require_set(params.V_0, "V_0");
  utils::require_set(params.K_0, "K_0");
  utils::require_set(params.Kprime_0, "Kprime_0");
  // Set defaults for missing values
  utils::fallback_to_default(params.E_0, 0.0);
  utils::fallback_to_default(params.P_0, 0.0);
  utils::fallback_to_default(params.G_0, std::nan(""));
  utils::fallback_to_default(params.Gprime_0, std::nan(""));
  // Check values are reasonable
  utils::check_in_range(params.V_0, 1.0e-7, 1.0e-3, "V_0");
  utils::check_in_range(params.K_0, 1.0e9, 1.0e13, "K_0");
  utils::check_in_range(params.Kprime_0, -5.0, 10.0, "Kprime_0");
}

// Compute P(V) - P as function to root find
double Vinet::vinet_gsl_wrapper(double x, void* p) {
  // Prefer C++ cast over C-style
  auto* vinet_params = static_cast<const gsl_params::SolverParams_P*>(p); // <const gsl_params::SolverParams_P*
  return compute_vinet(x / *(vinet_params->params.V_0), vinet_params->params)
    - vinet_params->pressure;
}
// For lambda can do:
// auto lambda = [&params, pressure](double x) {
//   return compute_vinet(x / params.V_0, params) - pressure;
// };

double Vinet::compute_vinet(
  double compression,
  const types::MineralParams& params
) {
  // Compute x^1/3 separately
  double eta = (3.0 / 2.0) * ((*params.Kprime_0) - 1.0);
  double x_cbrt = std::cbrt(compression);
  double mu = 1.0 - x_cbrt;
  return
    3.0 * (*params.K_0)
    * (1.0 / (x_cbrt * x_cbrt))
    * mu
    * std::exp(eta * mu)
    + (*params.P_0);
}

double Vinet::compute_volume(
  double pressure,
  double temperature [[maybe_unused]],
  const types::MineralParams& params
) const {
  // Find root in [P(V) - P] to find V that fits P
  // Set up GSL function params
  gsl_params::SolverParams_P vinet_params{params, pressure};
  // Set a, b limits
  double x_lo = 0.1 * (*params.V_0);
  double x_hi = 1.5 * (*params.V_0);
  double volume_root = optim::roots::brent(
    &vinet_gsl_wrapper,
    vinet_params,
    x_lo,
    x_hi
  );
  return volume_root;
}

double Vinet::compute_pressure(
  double temperature [[maybe_unused]],
  double volume,
  const types::MineralParams& params
) const {
  return compute_vinet(volume/(*params.V_0), params);
}

double Vinet::compute_isothermal_bulk_modulus_reuss(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume,
  const types::MineralParams& params
) const {
  double eta = (3.0/2.0) * ((*params.Kprime_0) - 1.0);
  double x_cbrt = std::cbrt(volume / (*params.V_0));
  double mu = 1.0 - x_cbrt;
  return
    ((*params.K_0) * (1.0 / (x_cbrt * x_cbrt)))
    * (1 + ((eta * x_cbrt + 1.0) * mu))
    * std::exp(eta * mu);
}

double Vinet::compute_isentropic_bulk_modulus_reuss(
  double pressure,
  double temperature,
  double volume,
  const types::MineralParams& params
) const {
  // K_T, K_S the same here
  return compute_isothermal_bulk_modulus_reuss(
    pressure, temperature, volume, params);
}

double Vinet::compute_gibbs_free_energy(
  double pressure,
  double temperature,
  double volume,
  const types::MineralParams& params
) const {
  // G = int VdP = [PV] - int PdV = E + PV
  return compute_molar_internal_energy(
    pressure, temperature, volume, params
    ) + volume * pressure;
}

double Vinet::compute_molar_internal_energy(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume,
  const types::MineralParams& params
) const {
  double eta = 1.5 * ((*params.Kprime_0) - 1.0);
  double x_cbrt = std::cbrt(volume / (*params.V_0));
  double mu = 1.0 - x_cbrt;
  double eta_mu = eta * mu;
  double intPdV = 9.0 * (*params.V_0) * (*params.K_0)
    / (eta * eta)
    * ((1.0 - eta_mu) * std::exp(eta_mu) - 1.0);
  return -intPdV + (*params.E_0);
}

// Meaningless parameters for isothermal EOS
// Returning defaults
double Vinet::compute_grueneisen_parameter(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double Vinet::compute_shear_modulus(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double Vinet::compute_molar_heat_capacity_v(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 1.0e99;
}

double Vinet::compute_molar_heat_capacity_p(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 1.0e99;
}

double Vinet::compute_thermal_expansivity(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double Vinet::compute_entropy(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

} // namespace burnman::eos
