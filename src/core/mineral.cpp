/*
  TODO: Copyright Notice!
*/
#include <cmath>
#include <stdexcept>
#include "burnman/core/mineral.hpp"
#include "burnman/eos/property_modifiers.hpp"
#include "burnman/utils/make_eos.hpp"
#include "burnman/utils/constants.hpp"

void Mineral::set_method(std::unique_ptr<EquationOfState> new_method) {
  // Here new_method is a unique pointer to an EOS class
  eos_method = std::move(new_method);
  // Set the params.equation_of_state to Custom
  params.equation_of_state = EOSType::Custom;
  // Clear material properties cache
  reset();
  // Validate parameters
  eos_method->validate_parameters(params);
}

void Mineral::set_method(EOSType new_method) {
  // Check if EOSType is Auto and update from params
  if (new_method == EOSType::Auto) {
    if (params.equation_of_state.has_value()) {
      new_method = *params.equation_of_state;
    } else {
      throw std::invalid_argument("No EOS set in parameters!");
    }
  }
  // Here the new_method is a predefined enum
  eos_method = make_eos(new_method);
  // Set the params.equation_of_state to new type
  params.equation_of_state = new_method;
  // Clear materials properties cache
  reset();
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
  property_modifier_excesses = excesses::Excesses();
  // Loop through param vector and call overloaded excess func
  for (const auto& excess_param : property_modifier_params) {
    std::visit([&](auto&& p) {
        property_modifier_excesses += excesses::compute_excesses(
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
  return (eos_method->compute_thermal_expansivity(
      get_pressure(),
      get_temperature(),
      get_molar_volume_unmodified(),
      params)
      + property_modifier_excesses.d2GdPdT)
    / get_molar_volume();
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

// TODO - formula

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
  double thermal_expansivity = get_thermal_expansivity();
  return
    get_molar_heat_capacity_p()
    - get_molar_volume()
    * get_temperature()
    * thermal_expansivity * thermal_expansivity
    * get_isothermal_bulk_modulus_reuss();
}

double Mineral::compute_isentropic_thermal_gradient() const {
  return
    (get_molar_volume() * get_temperature() * get_thermal_expansivity())
    / get_molar_heat_capacity_p();
}


// TODO:
//       implement set_method so know how to call EOS funcs

// Mineral parameter properties - util/eos MineralParams
// Property modifier params - util/eos ExcessParams::(Struct) --> use std::visit with overloaded functions to compute correct excess function
// Property modifiers - util/eos Excesses - has overloaed += operator to add data when computing loop

// Pass properties as value (unless string, vector, etc.)
// Pass params as const reference