/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/eos/components/debye.hpp"
#include <cmath>
#include <gsl/gsl_sf_debye.h>
#include <gsl/gsl_integration.h>
#include "burnman/utils/constants.hpp"

namespace burnman::eos {

double debye::debye_fn_integrand(double xi, void*) {
  return (xi * xi * xi) / std::expm1(xi);
}

double debye::debye_fn_quad(double x) {
  gsl_integration_workspace* gsl_w =
    gsl_integration_workspace_alloc(1000);
  double result, error;
  gsl_function integrand;
  integrand.function = &debye_fn_integrand;
  gsl_integration_qags(
    &integrand,
    0, x, // limits
    1.49e-8, 1.49e-8, //eps from scipy TODO: move to constants
    1000, // limit
    gsl_w, &result, &error);
  // error unused here
  // result hold integral
  // Free allocated memory
  gsl_integration_workspace_free(gsl_w);
  return 3.0 * result / (x*x*x);
}

// For errors do:
// gsl_sf_result result;
// int status = gsl_sf_debye_3_e(x, &result);
// for a string do: gsl_strerror(status);
// result.val, result.err
double debye::debye_fn_cheb(double x) {
  return gsl_sf_debye_3(x);
}

double debye::compute_thermal_energy_impl(
  double temperature,
  double debye_temperature,
  double napfu) {
  if (temperature <= constants::precision::double_eps) {
    return 0.0;
  }
  // E_th
  return 3.0 * napfu
    * constants::physics::gas_constant
    * temperature
    * debye_fn_cheb(debye_temperature / temperature);
}

double debye::compute_molar_heat_capacity_v_impl(
  double temperature,
  double debye_temperature,
  double napfu) {
  if (temperature < constants::precision::double_eps) {
    return 0.0;
  }
  double x = debye_temperature / temperature;
  // C_v
  return 3.0 * napfu
    * constants::physics::gas_constant
    * (4.0 * debye_fn_cheb(x) - 3.0 * x / std::expm1(x));
}

double debye::compute_helmholtz_free_energy_impl(
  double temperature,
  double debye_temperature,
  double napfu) {
  if (temperature <= constants::precision::double_eps) {
    return 0.0;
  }
  double x = debye_temperature / temperature;
  // F
  return napfu * constants::physics::gas_constant
    * temperature
    * (3.0 * std::log1p(-std::exp(-x)) - debye_fn_cheb(x));
}

double debye::compute_entropy_impl(
  double temperature,
  double debye_temperature,
  double napfu) {
  if (temperature <= constants::precision::double_eps) {
    return 0.0;
  }
  double x = debye_temperature / temperature;
  // S
  return napfu * constants::physics::gas_constant
    * (4.0 * debye_fn_cheb(x) - 3.0 * std::log1p(-std::exp(-x))
    );
}

double debye::compute_dmolar_heat_capacity_v_dT_impl(
  double temperature,
  double debye_temperature,
  double napfu) {
  if (temperature <= constants::precision::double_eps) {
    return 0.0;
  }
  double x = debye_temperature / temperature;
  double Cv_over_T = compute_molar_heat_capacity_v_impl(
    temperature,
    debye_temperature,
    napfu)
    / temperature;
  double E_over_Tsqr = compute_thermal_energy_impl(
    temperature,
    debye_temperature,
    napfu)
    / (temperature * temperature);
  // dCvdT
  return 3.0 * Cv_over_T
    + (Cv_over_T - 4.0 * E_over_Tsqr)
    * x / (1.0 - std::exp(-x));
}

double debye::compute_thermal_energy(
  double temperature,
  double debye_temperature,
  int napfu) {
  return compute_thermal_energy_impl(
    temperature,
    debye_temperature,
    static_cast<double>(napfu)
  );
}

double debye::compute_thermal_energy(
  double temperature,
  double debye_temperature,
  types::ExplicitDouble f) {
  return compute_thermal_energy_impl(
    temperature,
    debye_temperature,
    f.value
  );
}

double debye::compute_molar_heat_capacity_v(
  double temperature,
  double debye_temperature,
  int napfu) {
  return compute_molar_heat_capacity_v_impl(
    temperature,
    debye_temperature,
    static_cast<double>(napfu)
  );
}

double debye::compute_molar_heat_capacity_v(
  double temperature,
  double debye_temperature,
  types::ExplicitDouble f) {
  return compute_molar_heat_capacity_v_impl(
    temperature,
    debye_temperature,
    f.value
  );
}

double debye::compute_helmholtz_free_energy(
  double temperature,
  double debye_temperature,
  int napfu) {
  return compute_helmholtz_free_energy_impl(
    temperature,
    debye_temperature,
    static_cast<double>(napfu)
  );
}

double debye::compute_helmholtz_free_energy(
  double temperature,
  double debye_temperature,
  types::ExplicitDouble f) {
  return compute_helmholtz_free_energy_impl(
    temperature,
    debye_temperature,
    f.value
  );
}

double debye::compute_entropy(
  double temperature,
  double debye_temperature,
  int napfu) {
  return compute_entropy_impl(
    temperature,
    debye_temperature,
    static_cast<double>(napfu)
  );
}

double debye::compute_entropy(
  double temperature,
  double debye_temperature,
  types::ExplicitDouble f) {
  return compute_entropy_impl(
    temperature,
    debye_temperature,
    f.value
  );
}

double debye::compute_dmolar_heat_capacity_v_dT(
  double temperature,
  double debye_temperature,
  int napfu) {
  return compute_dmolar_heat_capacity_v_dT_impl(
    temperature,
    debye_temperature,
    static_cast<double>(napfu)
  );
}

double debye::compute_dmolar_heat_capacity_v_dT(
  double temperature,
  double debye_temperature,
  types::ExplicitDouble f) {
  return compute_dmolar_heat_capacity_v_dT_impl(
    temperature,
    debye_temperature,
    f.value
  );
}

} // namespace burnman::eos
