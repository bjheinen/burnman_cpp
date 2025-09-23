/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/core/equation_of_state.hpp"
#include <typeinfo>
#include "burnman/utils/exceptions.hpp"
#include "burnman/utils/warnings.hpp"

namespace burnman {

// Helper functions to throw errors
[[noreturn]] void EquationOfState::throw_not_implemented_error(const std::string& method) const {
  throw exceptions::NotImplementedError(get_class_name(), method);
}

std::string EquationOfState::get_class_name() const {
  return typeid(*this).name();
}

// Implemented EOS functions
double EquationOfState::compute_density(
  double volume,
  const types::MineralParams& params
) const {
  return *params.molar_mass / volume;
}

// Default compute function implementations (override in derived classes)

void EquationOfState::validate_parameters(types::MineralParams& params [[maybe_unused]]){
  // Pass with warning as default
  utils::warn("No parameter validation implemented for EOS: " + get_class_name());
}

double EquationOfState::compute_volume(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_pressure(
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_grueneisen_parameter(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_isothermal_bulk_modulus_reuss(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_isentropic_bulk_modulus_reuss(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_shear_modulus(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_molar_heat_capacity_v(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_molar_heat_capacity_p(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_thermal_expansivity(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_gibbs_free_energy(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_helmholtz_free_energy(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_entropy(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_enthalpy(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_molar_internal_energy(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  double volume [[maybe_unused]],
  const types::MineralParams& params [[maybe_unused]]
) const {
  throw_not_implemented_error(__func__);
}

} // namespace burnman
