/*
  TODO: Copyright Notice!
*/
#include <typeinfo>
#include "burnman/core/material.hpp"
#include "burnman/utils/exceptions.hpp"

[[noreturn]] void Material::throw_not_implemented_error(const std::string& method) const {
  throw NotImplementedError(get_class_name(), method);
}

std::string Material::get_class_name() const {
  return typeid(*this).name();
}

void Material::reset() {
  molar_internal_energy.reset();
  molar_gibbs.reset();
  molar_helmholtz.reset();
  molar_mass.reset();
  molar_volume.reset();
  molar_volume_unmodified.reset();
  density.reset();
  molar_entropy.reset();
  molar_enthalpy.reset();
  isothermal_bulk_modulus_reuss.reset();
  isentropic_bulk_modulus_reuss.reset();
  isothermal_compressibility_reuss.reset();
  isentropic_compressibility_reuss.reset();
  shear_modulus.reset();
  p_wave_velocity.reset();
  bulk_sound_velocity.reset();
  shear_wave_velocity.reset();
  grueneisen_parameter.reset();
  thermal_expansivity.reset();
  molar_heat_capacity_v.reset();
  molar_heat_capacity_p.reset();
  isentropic_thermal_gradient.reset();
}

void Material::set_name(std::string new_name) {
  name = new_name;
}

std::string Material::get_name() const {
  if (name.has_value()) {
    return *name;
  } else {
    return get_class_name();
  }
}

void Material::set_method(EOSType new_method) {
  throw_not_implemented_error(__func__);
}

void Material::set_method(std::unique_ptr<EquationOfState> new_method) {
  throw_not_implemented_error(__func__);
}

void Material::set_state(double new_pressure, double new_temperature) {
  reset();
  pressure = new_pressure;
  temperature = new_temperature;
}

double Material::get_pressure() const {
  return *pressure;
}

double Material::get_temperature() const {
  return *temperature;
}

// TODO: maybe factor out caching logic with function template?

double Material::get_molar_internal_energy() const {
  if (!molar_internal_energy.has_value()) {
    molar_internal_energy = compute_molar_internal_energy();
  }
  return *molar_internal_energy;
}

double Material::get_molar_gibbs() const {
  if (!molar_gibbs.has_value()) {
    molar_gibbs = compute_molar_gibbs();
  }
  return *molar_gibbs;
}

double Material::get_molar_helmholtz() const {
  if (!molar_helmholtz.has_value()) {
    molar_helmholtz = compute_molar_helmholtz();
  }
  return *molar_helmholtz;
}

double Material::get_molar_mass() const {
  if (!molar_mass.has_value()) {
    molar_mass = compute_molar_mass();
  }
  return *molar_mass;
}

double Material::get_molar_volume() const {
  if (!molar_volume.has_value()) {
    molar_volume = compute_molar_volume();
  }
  return *molar_volume;
}

double Material::get_molar_volume_unmodified() const {
  if (!molar_volume_unmodified.has_value()) {
    molar_volume_unmodified = compute_molar_volume_unmodified();
  }
  return *molar_volume_unmodified;
}

double Material::get_density() const {
  if (!density.has_value()) {
    density = compute_density();
  }
  return *density;
}

double Material::get_molar_entropy() const {
  if (!molar_entropy.has_value()) {
    molar_entropy = compute_molar_entropy();
  }
  return *molar_entropy;
}

double Material::get_molar_enthalpy() const {
  if (!molar_enthalpy.has_value()) {
    molar_enthalpy = compute_molar_enthalpy();
  }
  return *molar_enthalpy;
}

double Material::get_isothermal_bulk_modulus_reuss() const {
  if (!isothermal_bulk_modulus_reuss.has_value()) {
    isothermal_bulk_modulus_reuss = compute_isothermal_bulk_modulus_reuss();
  }
  return *isothermal_bulk_modulus_reuss;
}

double Material::get_isentropic_bulk_modulus_reuss() const {
  if (!isentropic_bulk_modulus_reuss.has_value()) {
    isentropic_bulk_modulus_reuss = compute_isentropic_bulk_modulus_reuss();
  }
  return *isentropic_bulk_modulus_reuss;
}

double Material::get_isothermal_compressibility_reuss() const {
  if (!isothermal_compressibility_reuss.has_value()) {
    isothermal_compressibility_reuss = compute_isothermal_compressibility_reuss();
  }
  return *isothermal_compressibility_reuss;
}

double Material::get_isentropic_compressibility_reuss() const {
  if (!isentropic_compressibility_reuss.has_value()) {
    isentropic_compressibility_reuss = compute_isentropic_compressibility_reuss();
  }
  return *isentropic_compressibility_reuss;
}

double Material::get_shear_modulus() const {
  if (!shear_modulus.has_value()) {
    shear_modulus = compute_shear_modulus();
  }
  return *shear_modulus;
}

double Material::get_p_wave_velocity() const {
  if (!p_wave_velocity.has_value()) {
    p_wave_velocity = compute_p_wave_velocity();
  }
  return *p_wave_velocity;
}

double Material::get_bulk_sound_velocity() const {
  if (!bulk_sound_velocity.has_value()) {
    bulk_sound_velocity = compute_bulk_sound_velocity();
  }
  return *bulk_sound_velocity;
}

double Material::get_shear_wave_velocity() const {
  if (!shear_wave_velocity.has_value()) {
    shear_wave_velocity = compute_shear_wave_velocity();
  }
  return *shear_wave_velocity;
}

double Material::get_grueneisen_parameter() const {
  if (!grueneisen_parameter.has_value()) {
    grueneisen_parameter = compute_grueneisen_parameter();
  }
  return *grueneisen_parameter;
}

double Material::get_thermal_expansivity() const {
  if (!thermal_expansivity.has_value()) {
    thermal_expansivity = compute_thermal_expansivity();
  }
  return *thermal_expansivity;
}

double Material::get_molar_heat_capacity_v() const {
  if (!molar_heat_capacity_v.has_value()) {
    molar_heat_capacity_v = compute_molar_heat_capacity_v();
  }
  return *molar_heat_capacity_v;
}

double Material::get_molar_heat_capacity_p() const {
  if (!molar_heat_capacity_p.has_value()) {
    molar_heat_capacity_p = compute_molar_heat_capacity_p();
  }
  return *molar_heat_capacity_p;
}

double Material::get_isentropic_thermal_gradient() const {
  if (!isentropic_thermal_gradient.has_value()) {
    isentropic_thermal_gradient = compute_isentropic_thermal_gradient();
  }
  return *isentropic_thermal_gradient;
}


// Default compute function implementations (override in derived classes)
double Material::compute_molar_internal_energy() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_molar_gibbs() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_molar_helmholtz() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_molar_mass() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_molar_volume() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_molar_volume_unmodified() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_density() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_molar_entropy() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_molar_enthalpy() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_isothermal_bulk_modulus_reuss() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_isentropic_bulk_modulus_reuss() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_isothermal_compressibility_reuss() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_isentropic_compressibility_reuss() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_shear_modulus() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_p_wave_velocity() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_bulk_sound_velocity() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_shear_wave_velocity() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_grueneisen_parameter() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_thermal_expansivity() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_molar_heat_capacity_v() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_molar_heat_capacity_p() const {
  throw_not_implemented_error(__func__);
}

double Material::compute_isentropic_thermal_gradient() const {
  throw_not_implemented_error(__func__);
}
