/*
  TODO: Copyright Notice!
*/
#include <cmath>
#include <stdexcept>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "burnman/eos/vinet.hpp"

// TODO: might need to access params.x with *

bool Vinet::validate_parameters(MineralParams& params) {

  // TODO: check if we can just do if (!params.V_0)

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
  if (*params.V_0 < 1.0e-7 || *params.V_0 > 1.0e-3) {
    ; //warning warnings.warn("Unusual value for V_0", stacklevel=2)
  }
  if (*params.K_0 < 1.0e9 || *params.K_0 > 1.0e13) {
    ; // warning warnings.warn("Unusual value for K_0", stacklevel=2)
  }
  if (*params.Kprime_0 < -5.0 || *params.Kprime_0 > 10.0) {
    ; // warning warnings.warn("Unusual value for Kprime_0", stacklevel=2)
  }

  return 1; // maybe make this void

}

// Compute P(V) - P as function to root find
double Vinet::vinet_gsl_wrapper(double x, void* p) {
  // Prefer C++ cast over C-style
  auto* vinet_params = static_cast<const ParamsGSL::SolverParams_P*>(p); // <const ParamsGSL::SolverParams_P*
  return compute_vinet(x / *(vinet_params->params.V_0), vinet_params->params)
    - vinet_params->pressure;
}
// For lambda can do:
// auto lambda = [&params, pressure](double x) {
//   return compute_vinet(x / params.V_0, params) - pressure;
// };

double Vinet::compute_vinet(
  double compression,
  const MineralParams& params
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
  const MineralParams& params
) const {
  // Find root in [P(V) - P] to find V that fits P
  // Set up GSL function object and params
  gsl_function vinet_functor;
  ParamsGSL::SolverParams_P vinet_params{params, pressure};
  vinet_functor.function = &vinet_gsl_wrapper;
  vinet_functor.params = &vinet_params;
  // Define and allocate the solver
  const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
  gsl_root_fsolver* solver = gsl_root_fsolver_alloc(T);
  // Set a, b limits
  double x_lo = 0.1 * (*params.V_0);
  double x_hi = 1.5 * (*params.V_0);
  // Set the GSL solver
  gsl_root_fsolver_set(solver, &vinet_functor, x_lo, x_hi);
  // Iterate the solver
  int iter = 0;
  int maxiter = 100;
  int status;
  double root;
  double atol = 2.0e-12;
  double rtol = 4.0 * std::numeric_limits<double>::epsilon();
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate(solver);
      root = gsl_root_fsolver_root(solver);
      x_lo = gsl_root_fsolver_x_lower(solver);
      x_hi = gsl_root_fsolver_x_upper(solver);
      status = gsl_root_test_interval(x_lo, x_hi,
                                      atol, rtol);
      // Leaving in print from GSL example to debug
      // TODO: remove prints when working/tests set up.
      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
              iter, x_lo, x_hi,
              root,
              x_hi - x_lo);
    } while (status == GSL_CONTINUE && iter < maxiter);

  // Get final root (can delet root getting above when tested)
  root = gsl_root_fsolver_root(solver);
  // Clean up
  gsl_root_fsolver_free(solver);

  if (status == GSL_SUCCESS) {
      return root;
  } else {
      return -1; // throw an error? "Error: Root-finding did not converge!"
  }
}

double Vinet::compute_pressure(
  double temperature [[maybe_unused]],
  double volume,
  const MineralParams& params
) const {
  return compute_vinet(volume/(*params.V_0), params);
}

double Vinet::compute_isothermal_bulk_modulus_reuss(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume,
  const MineralParams& params
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
  const MineralParams& params
) const {
  // K_T, K_S the same here
  return compute_isothermal_bulk_modulus_reuss(
    pressure, temperature, volume, params);
}

double Vinet::compute_gibbs_free_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
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
  const MineralParams& params
) const {
  double eta = (3.0/2.0) * ((*params.Kprime_0) - 1.0);
  double x_cbrt = std::cbrt(volume / (*params.V_0));
  double mu = 1.0 - x_cbrt;
  double intPdV = 9.0 * (*params.V_0) * (*params.K_0)
    / (eta * eta)
    * ((1.0 - eta * mu) * std::exp(eta * mu) - 1.0);
  return -intPdV + (*params.E_0);
}

// Meaningless parameters for isothermal EOS
// Returning defaults
double Vinet::compute_grueneisen_parameter(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double Vinet::compute_shear_modulus(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double Vinet::compute_molar_heat_capacity_v(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params [[maybe_unused]]
) const {
  return 1.0e99;
}

double Vinet::compute_molar_heat_capacity_p(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params [[maybe_unused]]
) const {
  return 1.0e99;
}

double Vinet::compute_thermal_expansivity(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}

double Vinet::compute_entropy(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const MineralParams& params [[maybe_unused]]
) const {
  return 0.0;
}
