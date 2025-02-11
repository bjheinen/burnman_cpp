/*
  TODO: Copyright Notice!
*/
#include <typeinfo>
#include "../../include/burnman/eos/equation_of_state.hpp"
#include "../../include/burnman/util/exceptions.hpp"

// Helper functions to throw errors
[[noreturn]] void EquationOfState::throw_not_implemented_error(const std::string& method) const {
  throw NotImplementedError(get_class_name(), method);
}

std::string EquationOfState::get_class_name() const {
  return typeid(*this).name();
}

// Implemented EOS functions
double EquationOfState::compute_density(
  double volume,
  const MineralParams& params
) const {
  return *params.molar_mass / volume;
}

// Default compute function implementations (override in derived classes)

bool EquationOfState::validate_parameters(MineralParams& params){
  return 1; // Pass quietly as default
}

double EquationOfState::compute_volume(
  double pressure,
  double temperature,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_pressure(
  double temperature,
  double volume,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_grueneisen_parameter(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_isothermal_bulk_modulus_reuss(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_isentropic_bulk_modulus_reuss(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_shear_modulus(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_molar_heat_capacity_v(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_molar_heat_capacity_p(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_thermal_expansivity(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_gibbs_free_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_helmholtz_free_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_entropy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_enthalpy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}

double EquationOfState::compute_molar_internal_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  throw_not_implemented_error(__func__);
}