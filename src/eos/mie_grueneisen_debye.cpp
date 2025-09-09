/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/eos/mie_grueneisen_debye.hpp"
#include <cmath>
#include <stdexcept>
#include "burnman/utils/constants.hpp"
#include "burnman/optim/roots.hpp"
#include "burnman/eos/models/debye.hpp"
#include "burnman/eos/birch_murnaghan.hpp"

bool MGD3::validate_parameters(MineralParams& params) {
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
  if (!params.molar_mass.has_value()) {
    throw std::invalid_argument("params object missing parameter: molar_mass");
  }
  if (!params.napfu.has_value()) {
    throw std::invalid_argument("params object missing parameter: napfu");
  }
  if (!params.debye_0.has_value()) {
    throw std::invalid_argument("params object missing parameter: debye_0");
  }
  if (!params.grueneisen_0.has_value()) {
    throw std::invalid_argument("params object missing parameter: grueneisen_0");
  }
  if (!params.q_0.has_value()) {
    throw std::invalid_argument("params object missing parameter: q_0");
  }
  // Set defaults for missing values
  if (!params.T_0.has_value()) {
    params.T_0 = 300.0;
  }
  if (!params.P_0.has_value()) {
    params.P_0 = 0.0;
  }
  if (!params.E_0.has_value()) {
    params.E_0 = 0.0;
  }
  if (!params.F_0.has_value()) {
    params.F_0 = 0.0;
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
  if (*params.T_0 < 0.0) {
    ; // warnings.warn("Unusual value for T_0", stacklevel=2)
  }
  if (*params.molar_mass < 0.001 || *params.molar_mass > 1.0) {
    ; //warning warnings.warn("Unusual value for molar_mass", stacklevel=2)
  }
  if (*params.napfu < 1 || *params.napfu > 1000) {
    ; // warning warnings.warn("Unusual value for napfu", stacklevel=2)
  }
  if (*params.debye_0 < 1.0 || *params.debye_0 > 10000.0) {
    ; // warning warnings.warn("Unusual value for debye_0", stacklevel=2)
  }
  if (*params.grueneisen_0 < 0.0 || *params.G_0 > 10.0) {
    ; // warnings.warn("Unusual value for grueneisen_0", stacklevel=2)
  }
  if (*params.q_0 < -10.0 || *params.q_0 > 10.0) {
    ; // warnings.warn("Unusual value for q_0", stacklevel=2)
  }
  return 1;
}

// Compute P(V) - P as function to root find
double MGD3::mgd_gsl_wrapper(double x, void* p) {
  auto* mgd_params = static_cast<const ParamsGSL::SolverParams_PT*>(p);
  return
    BM3::compute_birch_murnaghan(
      *mgd_params->params.V_0 / x,
      mgd_params->params)
    + compute_thermal_pressure(
      mgd_params->temperature,
      x,
      mgd_params->params)
    - compute_thermal_pressure(
      *mgd_params->params.T_0,
      x,
      mgd_params->params)
    - mgd_params->pressure;
}

double MGD3::compute_volume(
  double pressure,
  double temperature,
  const MineralParams& params
) const {
  // Find root in [P(V) - P] to find V that fits P
  // Set up GSL params struct - to pass to objective function
  ParamsGSL::SolverParams_PT mgd_params{
    params,
    pressure,
    temperature};
  // Set a, b limits
  double x_lo = 0.1 * (*params.V_0);
  double x_hi = 1.5 * (*params.V_0);
  double volume_root = roots::brent(
    &mgd_gsl_wrapper,
    mgd_params,
    x_lo,
    x_hi
  );
  return volume_root;
}

double MGD3::compute_pressure(
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return BM3::compute_birch_murnaghan(
    *params.V_0 / volume,
    params)
  + compute_thermal_pressure(temperature, volume, params)
  - compute_thermal_pressure(*params.T_0, volume, params);
}

double MGD3::compute_grueneisen_parameter(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume,
  const MineralParams& params
) const {
  return compute_mgd_grueneisen_parameter(
    *params.V_0 / volume,
    params);
}

double MGD3::compute_mgd_grueneisen_parameter(
  double x,
  const MineralParams& params
) {
  return *params.grueneisen_0 * std::pow(x, -*params.q_0);
}

double MGD3::compute_isothermal_bulk_modulus_reuss(
  double pressure [[maybe_unused]],
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return BM3::compute_bm_bulk_modulus(volume, params)
    + compute_thermal_bulk_modulus(temperature, volume, params)
    - compute_thermal_bulk_modulus(*params.T_0, volume, params);
}

double MGD3::compute_isentropic_bulk_modulus_reuss(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double K_T = compute_isothermal_bulk_modulus_reuss(
    pressure,
    temperature,
    volume,
    params
  );
  double alpha = compute_thermal_expansivity(
    pressure,
    temperature,
    volume,
    params
  );
  double gamma = compute_mgd_grueneisen_parameter(
    *params.V_0 / volume,
    params);
  return K_T * (1.0 + gamma * alpha * temperature);
  // alpha is: gamma * C_v / K / volume (refactor?)
}

// Second order (MGD2)
double MGD2::compute_shear_modulus(
  double pressure [[maybe_unused]],
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return BM2::compute_second_order_shear_modulus(volume, params)
    + compute_thermal_shear_modulus(temperature, volume, params)
    - compute_thermal_shear_modulus(*params.T_0, volume, params);
}

// Third order (MGD3)
double MGD3::compute_shear_modulus(
  double pressure [[maybe_unused]],
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return BM3::compute_third_order_shear_modulus(volume, params)
    + compute_thermal_shear_modulus(temperature, volume, params)
    - compute_thermal_shear_modulus(*params.T_0, volume, params);
}

double MGD3::compute_molar_heat_capacity_v(
  double pressure [[maybe_unused]],
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double debye_temperature = compute_debye_temperature(
    *params.V_0/volume,
    params);
  return debye::compute_molar_heat_capacity_v(
    temperature,
    debye_temperature,
    *params.napfu
  );
}

double MGD3::compute_molar_heat_capacity_p(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double C_v = compute_molar_heat_capacity_v(
    pressure,
    temperature,
    volume,
    params
  );
  double gamma = compute_mgd_grueneisen_parameter(
    *params.V_0 / volume,
    params
  );
  double alpha = compute_thermal_expansivity(
    pressure,
    temperature,
    volume,
    params
  );
  return C_v * (1.0 + gamma * alpha * temperature);
}

double MGD3::compute_thermal_expansivity(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double C_v = compute_molar_heat_capacity_v(
    pressure,
    temperature,
    volume,
    params
  );
  double gamma = compute_mgd_grueneisen_parameter(
    *params.V_0 / volume,
    params
  );
  double K = compute_isothermal_bulk_modulus_reuss(
    pressure,
    temperature,
    volume,
    params
  );
  return gamma * C_v / K / volume;
}

double MGD3::compute_gibbs_free_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return compute_helmholtz_free_energy(
    pressure,
    temperature,
    volume,
    params
  )
  + pressure * volume;
}

double MGD3::compute_entropy(
  double pressure [[maybe_unused]],
  double temperature,
  double volume,
  const MineralParams& params [[maybe_unused]]
) const {
  double debye_temperature = compute_debye_temperature(
    *params.V_0 / volume,
    params
  );
  return debye::compute_entropy(
    temperature,
    debye_temperature,
    *params.napfu);
}

double MGD3::compute_molar_internal_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return compute_helmholtz_free_energy(
    pressure,
    temperature,
    volume,
    params)
  + temperature
  * compute_entropy(
    pressure,
    temperature,
    volume,
    params);
}

double MGD3::compute_helmholtz_free_energy(
  double pressure [[maybe_unused]],
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double x = *params.V_0 / volume;
  double x_cbrt = std::cbrt(x);
  double f = 0.5 * ((x_cbrt*x_cbrt) - 1.0);
  // Eq. 28/29 in SLB2005
  double b_iikk = 9.0 * (*params.K_0);
  double b_iikkmm = 27.0 * (*params.K_0) * (*params.Kprime_0 - 4.0);
  constexpr double one_sixth = 1.0 / 6.0;
  double F_pressure =
    0.5 * b_iikk * f * f * (*params.V_0)
    + one_sixth * (*params.V_0) * b_iikkmm * f * f * f;
  double debye_temperature = compute_debye_temperature(x, params);
  double F_thermal =
    debye::compute_helmholtz_free_energy(
      temperature, debye_temperature, *params.napfu)
    - debye::compute_helmholtz_free_energy(
      *params.T_0, debye_temperature, *params.napfu);
  return *params.F_0 + F_pressure + F_thermal;
}

double MGD3::compute_enthalpy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return compute_helmholtz_free_energy(
    pressure,
    temperature,
    volume,
    params)
  + temperature
  * compute_entropy(
    pressure,
    temperature,
    volume,
    params)
  + pressure * volume;
}

double MGD3::compute_thermal_shear_modulus(
  double temperature,
  double volume,
  const MineralParams& params
) const {
  if (temperature > 1.0e-10) {
    double x = *params.V_0 / volume;
    double gamma = compute_mgd_grueneisen_parameter(x, params);
    double debye_temperature = compute_debye_temperature(x, params);
    return 0.6 * (
      compute_thermal_bulk_modulus(temperature, volume, params)
      - 6.0 * constants::physics::gas_constant * temperature * (*params.napfu)
      / volume
      * gamma
      * debye::debye_fn_cheb(debye_temperature / temperature)
      );
  } else {
    return 0.0;
  }
}

double MGD3::compute_debye_temperature(
  double x,
  const MineralParams& params
) {
  double gamma_diff = *params.grueneisen_0
    - compute_mgd_grueneisen_parameter(x,params);
  return *params.debye_0 * std::exp(gamma_diff / *params.q_0);
}

double MGD3::compute_thermal_pressure(
  double temperature,
  double volume,
  const MineralParams& params
) {
  double x = *params.V_0 / volume;
  double debye_temperature = compute_debye_temperature(x, params);
  double gamma = compute_mgd_grueneisen_parameter(x, params);
  double E_th = debye::compute_thermal_energy(
    temperature,
    debye_temperature,
    *params.napfu);
  return gamma * E_th / volume;
}

double MGD3::compute_thermal_bulk_modulus(
  double temperature,
  double volume,
  const MineralParams& params
) const {
  if (temperature > 1.0e-10) {
    double x = *params.V_0 / volume;
    double gamma = compute_mgd_grueneisen_parameter(x, params);
    double debye_temperature = compute_debye_temperature(x, params);
    double xi = debye_temperature / temperature;
    return 3.0 * (*params.napfu) * constants::physics::gas_constant
      * temperature / volume * gamma
      * (
        (1.0 - (*params.q_0) - 3.0 * gamma) * debye::debye_fn_cheb(xi)
        + 3.0 * gamma * xi / (std::exp(xi) - 1.0));
  } else {
    return 0.0;
  }
}
