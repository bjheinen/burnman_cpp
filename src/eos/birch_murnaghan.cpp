/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/eos/birch_murnaghan.hpp"
#include <cmath>
#include <stdexcept>
#include "burnman/eos/components/gsl_params.hpp"
#include "burnman/optim/roots/brent.hpp"

namespace burnman::eos {

bool BM3::validate_parameters(types::MineralParams& params) {

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

  // Set defaults for missing values
  if (!params.E_0.has_value()) {
    params.E_0 = 0.0;
  }
  if (!params.P_0.has_value()) {
    params.P_0 = 0.0;
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
  if (*params.V_0 < 1.0e-7 || *params.V_0 > 1.0e-3) {
    ; //warning warnings.warn("Unusual value for V_0", stacklevel=2)
  }
  if (*params.K_0 < 1.0e9 || *params.K_0 > 1.0e13) {
    ; // warning warnings.warn("Unusual value for K_0", stacklevel=2)
  }
  if (*params.Kprime_0 < 0.0 || *params.Kprime_0 > 20.0) {
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

// Compute P(V) - P as function to root find
double BM3::bm_gsl_wrapper(double x, void* p) {
  auto* bm_params = static_cast<const gsl_params::SolverParams_P*>(p);
  return compute_birch_murnaghan(*bm_params->params.V_0/x, bm_params->params)
    - bm_params->pressure;
}

double BM3::compute_birch_murnaghan(
  double inv_compression,
  const types::MineralParams& params
) {
  double cbrt_x = std::cbrt(inv_compression);
  double x23 = cbrt_x * cbrt_x;
  double x53 = inv_compression * x23;
  double x73 = inv_compression * inv_compression * cbrt_x;
  return 1.5 * (*params.K_0)
    * (x73 - x53)
    * (1 + (0.75 * (*params.Kprime_0 - 4) * (x23 - 1)));
}

double BM3::compute_volume(
  double pressure,
  double temperature [[maybe_unused]],
  const types::MineralParams& params
) const {
  // Find root in [P(V) - P] to find V that fits P
  // Set up GSL params struct - to pass to objective function
  gsl_params::SolverParams_P bm_params{params, pressure};
  // Set a, b limits
  double x_lo = 0.1 * (*params.V_0);
  double x_hi = 1.5 * (*params.V_0);
  double volume_root = optim::roots::brent(
    &bm_gsl_wrapper,
    bm_params,
    x_lo,
    x_hi
  );
  return volume_root;
}

double BM3::compute_pressure(
  double temperature [[maybe_unused]],
  double volume,
  const types::MineralParams& params
) const {
  // Prefer V_0/V (x^-1) for faster exp
  return compute_birch_murnaghan((*params.V_0)/volume, params);
}

double BM3::compute_bm_bulk_modulus(
  double volume,
  const types::MineralParams& params
) {
  double x = (*params.V_0) / volume;
  double x_cbrt = std::cbrt(x);
  double f = 0.5 * ((x_cbrt*x_cbrt) - 1.0);
  // (1 + 2f)^(5/2) = (1+2f)^2 * (1+2f)^(1/2)
  double t = 1 + 2*f;
  double K = t * t * std::sqrt(t) * (
    (*params.K_0)
    + (3.0 * (*params.K_0) * (*params.Kprime_0) - 5.0 * (*params.K_0)) * f
    + 13.5 * ((*params.K_0) * (*params.Kprime_0) - 4.0 * (*params.K_0)) * f * f
  );
  return K;
}

double BM3::compute_isothermal_bulk_modulus_reuss(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume,
  const types::MineralParams& params
) const {
  return compute_bm_bulk_modulus(volume, params);
}

double BM3::compute_isentropic_bulk_modulus_reuss(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume,
  const types::MineralParams& params
) const {
  return compute_bm_bulk_modulus(volume, params);
}

double BM3::compute_gibbs_free_energy(
  double pressure,
  double temperature,
  double volume,
  const types::MineralParams& params
) const {
  return volume * pressure
    + compute_molar_internal_energy(
        pressure, temperature, volume, params);
}

double BM3::compute_molar_internal_energy(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume,
  const types::MineralParams& params
) const {
  double x = std::cbrt(*params.V_0 / volume);
  x *= x; //(V_0/V)^(2/3)
  double t = x - 1;
  // 9/16 = 0.5625
  double fac = 0.5625 * (*params.V_0) * (*params.K_0);
  double intPdV = fac * ((t*t*t * (*params.Kprime_0)) + (t*t * (6 - 4*x)));
  return *params.E_0 + intPdV;
}

double BM2::compute_second_order_shear_modulus(
  double volume,
  const types::MineralParams& params
) {
  double x = *params.V_0 / volume;
  double cbrt_x = std::cbrt(x);
  double x23 = cbrt_x * cbrt_x;
  // double x53 = x * x23;
  double G = (*params.G_0)
    * x * x23
    * (
      1.0
      - 0.5
      * (x23 - 1.0)
      * (5.0 - 3.0 * (*params.Gprime_0) * (*params.K_0) / (*params.G_0))
    );
    return G;
}

double BM3::compute_third_order_shear_modulus(
  double volume,
  const types::MineralParams& params
) {
  double x = *params.V_0 / volume;
  double x_cbrt = std::cbrt(x);
  double f = 0.5 * ((x_cbrt*x_cbrt) - 1.0);
  // (1 + 2f)^(5/2) = (1+2f)^2 * (1+2f)^(1/2)
  double t = 1 + 2*f;
  double G = t * t * std::sqrt(t)
    * (
      (*params.G_0)
      + (3.0 * (*params.K_0) * (*params.Gprime_0) - 5.0 * (*params.G_0))
      * f
      + (
        6.0 * (*params.K_0) * (*params.Gprime_0)
        - 24.0 * (*params.K_0)
        - 14.0 * (*params.G_0)
        + 4.5 * (*params.K_0) * (*params.Kprime_0))
      * f * f
    );
  return G;
}

// Second order (BM2)
double BM2::compute_shear_modulus(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume,
  const types::MineralParams& params
) const {
  return compute_second_order_shear_modulus(volume, params);
}

// Third order (BM3)
double BM3::compute_shear_modulus(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume,
  const types::MineralParams& params
) const {
  return compute_third_order_shear_modulus(volume, params);
}

// Meaningless parameters for isothermal EOS
// Returning defaults
double BM3::compute_grueneisen_parameter(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double BM3::compute_molar_heat_capacity_v(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 1.0e99;
}

double BM3::compute_molar_heat_capacity_p(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 1.0e99;
}

double BM3::compute_thermal_expansivity(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double BM3::compute_entropy(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

} // namespace burnman::eos
