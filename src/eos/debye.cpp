/*
  TODO: Copyright Notice!
*/
#include <cmath>
//#include <stdexcept>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_debye.h>
#include <gsl/gsl_integration.h>
#include "burnman/eos/debye.hpp"
#include "burnman/utils/constants.hpp"

double Debye::debye_fn_integrand(double xi, void*) {
  return (xi * xi * xi) / (std::exp(xi) - 1.0);
}

double Debye::debye_fn_quad(double x) {
  gsl_integration_workspace* gsl_w =
    gsl_integration_workspace_alloc(1000);
  double result, error;
  gsl_function integrand;
  integrand.function = &debye_fn_integrand;
  gsl_integration_qags(
    &integrand,
    0, x, // limits
    1.49e-8, 1.49e-8, //eps from scipy
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
double Debye::debye_fn_cheb(double x) {
  return gsl_sf_debye_3(x);
}

double Debye::compute_thermal_energy(
  double temperature,
  double debye_temperature,
  int napfu) {
  if (temperature <= constants::precision::double_eps) {
    return 0.0;
  }
  // E_th
  return 3.0 * napfu
    * constants::physics::gas_constant
    * temperature
    * debye_fn_cheb(debye_temperature / temperature);
}


double Debye::compute_molar_heat_capacity_v(
  double temperature,
  double debye_temperature,
  int napfu) {
  if (temperature < constants::precision::double_eps) {
    return 0.0;
  }
  double x = debye_temperature / temperature;
  // C_v
  return 3.0 * napfu
    * constants::physics::gas_constant
    * (4.0 * debye_fn_cheb(x) - 3.0 * x / (std::exp(x) - 1.0));
}

double Debye::compute_helmholtz_free_energy(
  double temperature,
  double debye_temperature,
  int napfu) {
  if (temperature <= constants::precision::double_eps) {
    return 0.0;
  }
  double x = debye_temperature / temperature;
  // F
  return napfu * constants::physics::gas_constant
    * temperature
    * (3.0 * std::log1p(-std::exp(-x)) - debye_fn_cheb(x));
}

double Debye::compute_entropy(
  double temperature,
  double debye_temperature,
  int napfu) {
  if (temperature <= constants::precision::double_eps) {
    return 0.0;
  }
  double x = debye_temperature / temperature;
  return napfu * constants::physics::gas_constant
    * (4.0 * debye_fn_cheb(x) - 3.0 * std::log1p(-std::exp(-x))
    );
}

double Debye::dmolar_heat_capacity_v_dT(
  double temperature,
  double debye_temperature,
  int napfu) {
  if (temperature <= constants::precision::double_eps) {
    return 0.0;
  }
  double x = debye_temperature / temperature;
  double Cv_over_T = compute_molar_heat_capacity_v(
    temperature,
    debye_temperature,
    napfu)
    / temperature;
  double E_over_Tsqr = compute_thermal_energy(
    temperature,
    debye_temperature,
    napfu)
    / (temperature * temperature);
  return 3.0 * Cv_over_T
    + (Cv_over_T - 4.0 * E_over_Tsqr)
    * x / (1.0 - std::exp(-x));
}

