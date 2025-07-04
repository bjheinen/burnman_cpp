/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
//#include <cmath>
#include "burnman/core/composite.hpp"
//#include "burnman/utils/constants.hpp"

void Composite::reset() {
  // Reset caches Material properties
  Material::reset();
  // Reset cached Composite properties
  volume_fractions.reset();
}

// Public getter functions with caching
Eigen::ArrayXd Composite::get_volume_fractions() const {
  if (!volume_fractions.has_value()) {
    volume_fractions = compute_volume_fractions();
  }
  return *volume_fractions;
}

// Material property overrides
double Composite::compute_molar_internal_energy() const {
 return (map_phases_to_array(&Material::get_molar_internal_energy)
   * molar_fractions).sum();
}

double Composite::compute_molar_gibbs() const {
  return (map_phases_to_array(&Material::get_molar_gibbs)
    * molar_fractions).sum();
}

double Composite::compute_molar_helmholtz() const {
  return (map_phases_to_array(&Material::get_molar_helmholtz)
    * molar_fractions).sum();
}

double Composite::compute_molar_mass() const {
  return (map_phases_to_array(&Material::get_molar_mass)
    * molar_fractions).sum();
}

double Composite::compute_molar_volume() const {
  return (map_phases_to_array(&Material::get_molar_volume)
    * molar_fractions).sum();
}

double Composite::compute_density() const {
  Eigen::ArrayXd densities = map_phases_to_array(&Material::get_density);
  Eigen::ArrayXd volumes = get_volume_fractions();
  return averaging_scheme->average_density(volumes, densities);
}

double Composite::compute_molar_entropy() const {
  return (map_phases_to_array(&Material::get_molar_entropy)
    * molar_fractions).sum();
}

double Composite::compute_molar_enthalpy() const {
  return (map_phases_to_array(&Material::get_molar_enthalpy)
    * molar_fractions).sum();
}

double Composite::compute_isothermal_bulk_modulus_reuss() const {
  Eigen::ArrayXd V_frac = get_volume_fractions();
  Eigen::ArrayXd K_ph = map_phases_to_array(&Material::get_isothermal_bulk_modulus_reuss);
  Eigen::ArrayXd G_ph = map_phases_to_array(&Material::get_shear_modulus);
  return averaging_scheme->average_bulk_moduli(V_frac, K_ph, G_ph);
}

double Composite::compute_isentropic_bulk_modulus_reuss() const {
  Eigen::ArrayXd V_frac = get_volume_fractions();
  Eigen::ArrayXd K_ph = map_phases_to_array(&Material::get_isentropic_bulk_modulus_reuss);
  Eigen::ArrayXd G_ph = map_phases_to_array(&Material::get_shear_modulus);
  return averaging_scheme->average_bulk_moduli(V_frac, K_ph, G_ph);
}

double Composite::compute_isothermal_compressibility_reuss() const {
  return 1.0 / get_isothermal_bulk_modulus_reuss();
}

double Composite::compute_isentropic_compressibility_reuss() const {
  return 1.0 / get_isentropic_bulk_modulus_reuss();
}

double Composite::compute_shear_modulus() const {
  Eigen::ArrayXd V_frac = get_volume_fractions();
  Eigen::ArrayXd K_ph = map_phases_to_array(&Material::get_isentropic_bulk_modulus_reuss);
  Eigen::ArrayXd G_ph = map_phases_to_array(&Material::get_shear_modulus);
  return averaging_scheme->average_shear_moduli(V_frac, K_ph, G_ph);
}

// TODO: think about redunancy/replication with `Solution::compute_p_wave...'
double Composite::compute_p_wave_velocity() const {
  constexpr double FOUR_THIRDS = 4.0 / 3.0;
  return std::sqrt(
    (get_isentropic_bulk_modulus_reuss() + FOUR_THIRDS * get_shear_modulus())
    / get_density()
  );
}

double Composite::compute_bulk_sound_velocity() const {
  return std::sqrt(
    get_isentropic_bulk_modulus_reuss() / get_density()
  );
}

double Composite::compute_shear_wave_velocity() const {
  return std::sqrt(
    get_shear_modulus() / get_density()
  );
}

double Composite::compute_grueneisen_parameter() const {
  return get_thermal_expansivity()
    * get_isothermal_bulk_modulus_reuss()
    * get_molar_volume()
    / get_molar_heat_capacity_v();
}

double Composite::compute_thermal_expansivity() const {
  Eigen::ArrayXd volumes = get_volume_fractions();
  Eigen::ArrayXd alphas = map_phases_to_array(&Material::get_thermal_expansivity);
  return averaging_scheme->average_thermal_expansivity(volumes, alphas);
}

double Composite::compute_molar_heat_capacity_v() const {
  Eigen::ArrayXd c_v = map_phases_to_array(&Material::get_molar_heat_capacity_v);
  return averaging_scheme->average_heat_capacity_v(molar_fractions, c_v);
}

double Composite::compute_molar_heat_capacity_p() const {
  Eigen::ArrayXd c_p = map_phases_to_array(&Material::get_molar_heat_capacity_p);
  return averaging_scheme->average_heat_capacity_p(molar_fractions, c_p);
}

// Composite properties
Eigen::ArrayXd Composite::compute_volume_fractions() const {
  return map_phases_to_array(&Material::get_molar_volume) * molar_fractions;
}
