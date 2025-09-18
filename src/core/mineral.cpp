/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/core/mineral.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>
#include <variant>
#include "burnman/utils/constants.hpp"
#include "burnman/eos/components/property_modifiers.hpp"
#include "burnman/eos/make_eos.hpp"

namespace burnman {

std::string Mineral::get_name() const {
  if (has_custom_name()) {
    return Material::get_name();
  } else if (params.name.has_value()) {
    //set_name(*params.name); // Need to deal with constness to set...
    return *params.name;
  } else {
    return Material::get_name();
  }
}

void Mineral::set_property_modifier_params(
  eos::excesses::ExcessParamVector excess_params
) {
  property_modifier_params = excess_params;
}

eos::excesses::Excesses Mineral::get_property_modifiers() const {
  return property_modifier_excesses;
}

void Mineral::set_method(std::shared_ptr<EquationOfState> new_method) {
  // Here new_method is a shared pointer to an EOS class
  eos_method = std::move(new_method);
  // Set the params.equation_of_state to Custom
  params.equation_of_state = types::EOSType::Custom;
  // Clear material properties cache
  reset_cache();
  // Validate parameters
  eos_method->validate_parameters(params);
}

void Mineral::set_method(types::EOSType new_method) {
  // Check if EOSType is Auto and update from params
  if (new_method == types::EOSType::Auto) {
    if (params.equation_of_state.has_value()) {
      new_method = *params.equation_of_state;
    } else {
      throw std::invalid_argument("No EOS set in parameters!");
    }
  }
  // Here the new_method is a predefined enum
  eos_method = eos::make_eos(new_method);
  // Set the params.equation_of_state to new type
  params.equation_of_state = new_method;
  // Clear materials properties cache
  reset_cache();
  // Validate parameters
  eos_method->validate_parameters(params);
}

void Mineral::set_state(
  double new_pressure,
  double new_temperature
) {
  // Set P,T using Material
  Material::set_state(new_pressure, new_temperature);
  // Compute property modifiers
  compute_property_modifiers();
}

void Mineral::compute_property_modifiers() {
  // Reset modifiers to zero
  property_modifier_excesses = eos::excesses::Excesses();
  // Loop through param vector and call overloaded excess func
  for (const auto& excess_param : property_modifier_params) {
    std::visit([&](auto&& p) {
        property_modifier_excesses += eos::excesses::compute_excesses(
          get_pressure(),
          get_temperature(),
          p);
    }, excess_param);
  }
}

// EOS properties - in P,T form
double Mineral::compute_molar_gibbs() const {
  return eos_method->compute_gibbs_free_energy(
    get_pressure(),
    get_temperature(),
    get_molar_volume_unmodified(),
    params)
    + property_modifier_excesses.G;
}

double Mineral::compute_molar_volume_unmodified() const {
  return eos_method->compute_volume(
    get_pressure(), get_temperature(), params);
}

double Mineral::compute_molar_volume() const {
  return get_molar_volume_unmodified() + property_modifier_excesses.dGdP;
}

double Mineral::compute_molar_entropy() const {
  return eos_method->compute_entropy(
    get_pressure(),
    get_temperature(),
    get_molar_volume_unmodified(),
    params)
    - property_modifier_excesses.dGdT;
}

double Mineral::compute_isothermal_bulk_modulus_reuss() const {
  double K_T_orig = eos_method->compute_isothermal_bulk_modulus_reuss(
    get_pressure(),
    get_temperature(),
    get_molar_volume_unmodified(),
    params);

  return get_molar_volume()
    / ((get_molar_volume_unmodified() / K_T_orig)
      - property_modifier_excesses.d2GdP2);
}

double Mineral::compute_molar_heat_capacity_p() const {
  return eos_method->compute_molar_heat_capacity_p(
    get_pressure(),
    get_temperature(),
    get_molar_volume_unmodified(),
    params)
    - get_temperature() * property_modifier_excesses.d2GdT2;
}

double Mineral::compute_thermal_expansivity() const {
  return (
    eos_method->compute_thermal_expansivity(
      get_pressure(),
      get_temperature(),
      get_molar_volume_unmodified(),
      params
    ) * get_molar_volume_unmodified()
    + property_modifier_excesses.d2GdPdT
  ) / get_molar_volume();
}

double Mineral::compute_shear_modulus() const {
  double G = eos_method->compute_shear_modulus(
    get_pressure(),
    get_temperature(),
    get_molar_volume_unmodified(),
    params);
  // TODO: Warning if G is 0 *(< machine eps)
  return G;
}

// Mineral properties (non EOS)
types::FormulaMap Mineral::compute_formula() const {
  if (params.formula.has_value()) {
    return *params.formula;
  } else {
    // can't use get_class_name at the moment because it is private
    throw std::invalid_argument("No Mineral formula set!");
  }
}

double Mineral::compute_molar_mass() const {
  if (params.molar_mass) {
    return *params.molar_mass;
  } else {
    throw std::invalid_argument("No molar mass parameter!");
  }
}

double Mineral::compute_density() const {
  return get_molar_mass() / get_molar_volume();
}

double Mineral::compute_molar_internal_energy() const {
  return
    get_molar_gibbs()
    - get_pressure() * get_molar_volume()
    + get_temperature() * get_molar_entropy();
}

double Mineral::compute_molar_helmholtz() const {
  return get_molar_gibbs() - get_pressure() * get_molar_volume();
}

double Mineral::compute_molar_enthalpy() const {
  return get_molar_gibbs() + get_temperature() * get_molar_entropy();
}

double Mineral::compute_isentropic_bulk_modulus_reuss() const {
  if (get_temperature() < 1.0e-10) {
    return get_isothermal_bulk_modulus_reuss();
  } else {
    return
      get_isothermal_bulk_modulus_reuss() * get_molar_heat_capacity_p()
      / get_molar_heat_capacity_v();
  }
}

double Mineral::compute_isothermal_compressibility_reuss() const {
  return 1.0 / get_isothermal_bulk_modulus_reuss();
}

double Mineral::compute_isentropic_compressibility_reuss() const {
  return 1.0 / get_isentropic_bulk_modulus_reuss();
}

double Mineral::compute_p_wave_velocity() const {
  return std::sqrt(
    (get_isentropic_bulk_modulus_reuss() + 4.0 / 3.0 * get_shear_modulus())
    / get_density()
  );
}

double Mineral::compute_bulk_sound_velocity() const {
  return std::sqrt(get_isentropic_bulk_modulus_reuss() / get_density());
}

double Mineral::compute_shear_wave_velocity() const {
  return std::sqrt(get_shear_modulus() / get_density());
}

double Mineral::compute_grueneisen_parameter() const {
  constexpr double eps = constants::precision::double_eps;
  if (std::fabs(get_molar_heat_capacity_v()) > eps) {
    return
      get_thermal_expansivity()
      * get_isothermal_bulk_modulus_reuss()
      * get_molar_volume()
      / get_molar_heat_capacity_v();
  } else if (
    std::fabs(property_modifier_excesses.d2GdPdT) < eps &&
    std::fabs(property_modifier_excesses.d2GdP2) < eps &&
    std::fabs(property_modifier_excesses.dGdP) < eps &&
    std::fabs(property_modifier_excesses.d2GdT2) < eps
    ) {
      return eos_method->compute_grueneisen_parameter(
        get_pressure(),
        get_temperature(),
        get_molar_volume(),
        params);
    } else {
      throw std::logic_error(
        "Error! You are trying to calculate the grueneisen parameter\n"
        "at a temperature where the heat capacity is very low and\n"
        "where you have defined Gibbs property modifiers.");
    }
}

double Mineral::compute_molar_heat_capacity_v() const {
  double alpha = get_thermal_expansivity();
  return
    get_molar_heat_capacity_p()
    - get_molar_volume()
    * get_temperature()
    * alpha * alpha
    * get_isothermal_bulk_modulus_reuss();
}

double Mineral::compute_isentropic_thermal_gradient() const {
  return
    (get_molar_volume() * get_temperature() * get_thermal_expansivity())
    / get_molar_heat_capacity_p();
}

} // namespace burnman
