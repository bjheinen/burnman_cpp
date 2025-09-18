/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/core/solution.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>
#include "burnman/utils/constants.hpp"
#include "burnman/utils/chemistry_utils.hpp"
#include "burnman/utils/matrix_utils.hpp"
#include "burnman/tools/averaging/averaging_schemes.hpp"

namespace burnman {

void Solution::reset_cache() {
  // Reset cached Material & CompositeMaterial properties
  CompositeMaterial::reset_cache();
  // Reset cached Solution properties
  excess_gibbs.reset();
  excess_volume.reset();
  excess_entropy.reset();
  excess_enthalpy.reset();
  activities.reset();
  activity_coefficients.reset();
  excess_partial_gibbs.reset();
  excess_partial_volumes.reset();
  excess_partial_entropies.reset();
  // partial_gibbs reset by CompositeMaterial::reset_cache()
  partial_volumes.reset();
  partial_entropies.reset();
  gibbs_hessian.reset();
  entropy_hessian.reset();
  volume_hessian.reset();
}

// Public setters
void Solution::set_solution_model(std::shared_ptr<solution_models::SolutionModel> model) {
  // TODO: think about using clone semantics instead
  solution_model = std::move(model);
}

void Solution::set_composition(const Eigen::ArrayXd& composition_vector) {
  // Throw error if no solution model yet
  if (!solution_model) {
    throw std::runtime_error("Cannot set molar fractions: solution model not set!");
  }
  // Throw error if length not correct
  if (static_cast<int>(composition_vector.size()) != get_n_endmembers()) {
    throw std::runtime_error(
      "Composition vector length (" + std::to_string(composition_vector.size()) +
      ") does not match number of endmembers (" + std::to_string(get_n_endmembers()) + ").");
  }
  if (std::abs(composition_vector.sum() - 1.0) > constants::precision::abs_tolerance) {
    throw std::runtime_error("Sum of molar fractions not equal to 1.0!");
  }
  molar_fractions = composition_vector;
}

// Public setter overrides of Material
void Solution::set_method(std::shared_ptr<EquationOfState> new_method) {
  for (Mineral& embr : solution_model->endmembers) {
    embr.set_method(new_method);
  }
  // Clear properties cache
  reset_cache();
}

void Solution::set_method(types::EOSType new_method) {
  for (Mineral& embr : solution_model->endmembers) {
    embr.set_method(new_method);
  }
  // Clear properties cache
  reset_cache();
}

void Solution::set_state(
  double new_pressure,
  double new_temperature
) {
  reset_cache(); // TODO: check if reset needed here???
  // Set P,T using Material
  Material::set_state(new_pressure, new_temperature);
  // Set state of each endmember
  for (Mineral& embr : solution_model->endmembers) {
    embr.set_state(new_pressure, new_temperature);
  }
}

Eigen::ArrayXd Solution::get_molar_fractions() const {
  return this->molar_fractions;
}

Eigen::ArrayXXd Solution::get_endmember_occupancies() const {
  return this->solution_model->endmember_occupancies;
}

Eigen::ArrayXXd Solution::get_endmember_n_occupancies() const {
  return this->solution_model->endmember_n_occupancies;
}

std::vector<std::string> Solution::get_site_names() const {
  return this->solution_model->site_names;
}

// Public getter functions with caching
double Solution::get_excess_gibbs() const {
  if (!excess_gibbs.has_value()) {
    excess_gibbs = compute_excess_gibbs();
  }
  return *excess_gibbs;
}

double Solution::get_excess_volume() const {
  if (!excess_volume.has_value()) {
    excess_volume = compute_excess_volume();
  }
  return *excess_volume;
}

double Solution::get_excess_entropy() const {
  if (!excess_entropy.has_value()) {
    excess_entropy = compute_excess_entropy();
  }
  return *excess_entropy;
}

double Solution::get_excess_enthalpy() const {
  if (!excess_enthalpy.has_value()) {
    excess_enthalpy = compute_excess_enthalpy();
  }
  return *excess_enthalpy;
}

Eigen::ArrayXd Solution::get_activities() const {
  if (!activities.has_value()) {
    activities = compute_activities();
  }
  return *activities;
}

Eigen::ArrayXd Solution::get_activity_coefficients() const {
  if (!activity_coefficients.has_value()) {
    activity_coefficients = compute_activity_coefficients();
  }
  return *activity_coefficients;
}

Eigen::ArrayXd Solution::get_excess_partial_gibbs() const {
  if (!excess_partial_gibbs.has_value()) {
    excess_partial_gibbs = compute_excess_partial_gibbs();
  }
  return *excess_partial_gibbs;
}

Eigen::ArrayXd Solution::get_excess_partial_volumes() const {
  if (!excess_partial_volumes.has_value()) {
    excess_partial_volumes = compute_excess_partial_volumes();
  }
  return *excess_partial_volumes;
}

Eigen::ArrayXd Solution::get_excess_partial_entropies() const {
  if (!excess_partial_entropies.has_value()) {
    excess_partial_entropies = compute_excess_partial_entropies();
  }
  return *excess_partial_entropies;
}

Eigen::ArrayXd Solution::get_partial_volumes() const {
  if (!partial_volumes.has_value()) {
    partial_volumes = compute_partial_volumes();
  }
  return *partial_volumes;
}

Eigen::ArrayXd Solution::get_partial_entropies() const {
  if (!partial_entropies.has_value()) {
    partial_entropies = compute_partial_entropies();
  }
  return *partial_entropies;
}

Eigen::MatrixXd Solution::get_gibbs_hessian() const {
  if (!gibbs_hessian.has_value()) {
    gibbs_hessian = compute_gibbs_hessian();
  }
  return *gibbs_hessian;
}

Eigen::MatrixXd Solution::get_entropy_hessian() const {
  if (!entropy_hessian.has_value()) {
    entropy_hessian = compute_entropy_hessian();
  }
  return *entropy_hessian;
}

Eigen::MatrixXd Solution::get_volume_hessian() const {
  if (!volume_hessian.has_value()) {
    volume_hessian = compute_volume_hessian();
  }
  return *volume_hessian;
}

// Solution property computations
double Solution::compute_excess_gibbs() const {
  return solution_model->compute_excess_gibbs_free_energy(
    get_pressure(), get_temperature(), molar_fractions);
}

double Solution::compute_excess_volume() const {
  return solution_model->compute_excess_volume(
    get_pressure(), get_temperature(), molar_fractions);
}

double Solution::compute_excess_entropy() const {
  return solution_model->compute_excess_entropy(
    get_pressure(), get_temperature(), molar_fractions);
}

double Solution::compute_excess_enthalpy() const {
  return solution_model->compute_excess_enthalpy(
    get_pressure(), get_temperature(), molar_fractions);
}

Eigen::ArrayXd Solution::compute_activities() const {
  return solution_model->compute_activities(
    get_pressure(), get_temperature(), molar_fractions);
}

Eigen::ArrayXd Solution::compute_activity_coefficients() const {
  return solution_model->compute_activity_coefficients(
    get_pressure(), get_temperature(), molar_fractions);
}

Eigen::ArrayXd Solution::compute_excess_partial_gibbs() const {
  return solution_model->compute_excess_partial_gibbs_free_energies(
    get_pressure(), get_temperature(), molar_fractions);
}

Eigen::ArrayXd Solution::compute_excess_partial_volumes() const {
  return solution_model->compute_excess_partial_volumes(
    get_pressure(), get_temperature(), molar_fractions);
}

Eigen::ArrayXd Solution::compute_excess_partial_entropies() const {
  return solution_model->compute_excess_partial_entropies(
    get_pressure(), get_temperature(), molar_fractions);
}

Eigen::ArrayXd Solution::compute_partial_gibbs() const {
  return get_excess_partial_gibbs()
    + map_endmembers_to_array(&Mineral::get_molar_gibbs);
}

Eigen::ArrayXd Solution::compute_partial_volumes() const {
  return get_excess_partial_volumes()
    + map_endmembers_to_array(&Mineral::get_molar_volume);
}

Eigen::ArrayXd Solution::compute_partial_entropies() const {
  return get_excess_partial_entropies()
    + map_endmembers_to_array(&Mineral::get_molar_entropy);
}

Eigen::MatrixXd Solution::compute_gibbs_hessian() const {
  return solution_model->compute_gibbs_hessian(
    get_pressure(), get_temperature(), molar_fractions);
}

Eigen::MatrixXd Solution::compute_entropy_hessian() const {
  return solution_model->compute_entropy_hessian(
    get_pressure(), get_temperature(), molar_fractions);
}

Eigen::MatrixXd Solution::compute_volume_hessian() const {
  return solution_model->compute_volume_hessian(
    get_pressure(), get_temperature(), molar_fractions);
}

// Material property overrides
double Solution::compute_molar_internal_energy() const {
  return get_molar_helmholtz()
    + get_temperature() * get_molar_entropy();
}

double Solution::compute_molar_gibbs() const {
  //Eigen::ArrayXd em_gibbs = map_endmembers_to_array(&Mineral::get_molar_gibbs);
  //return (em_gibbs * molar_fractions).sum() + get_excess_gibbs();
  return (get_partial_gibbs() * molar_fractions).sum();
}

double Solution::compute_molar_helmholtz() const {
  return get_molar_gibbs()
    - get_pressure() * get_molar_volume();
}

double Solution::compute_molar_mass() const {
  Eigen::ArrayXd masses = map_endmembers_to_array(&Mineral::get_molar_mass);
  return (masses * molar_fractions).sum();
}

double Solution::compute_molar_volume() const {
  return (get_partial_volumes() * molar_fractions).sum();
}

double Solution::compute_density() const {
  return get_molar_mass() / get_molar_volume();
}

double Solution::compute_molar_entropy() const {
  return (get_partial_entropies() * molar_fractions).sum();
}

double Solution::compute_molar_enthalpy() const {
  Eigen::ArrayXd em_enthalpies = map_endmembers_to_array(&Mineral::get_molar_enthalpy);
  return (em_enthalpies * molar_fractions).sum() + get_excess_enthalpy();
}

double Solution::compute_isothermal_bulk_modulus_reuss() const {
  double V_over_KT = (
    map_endmembers_to_array(&Mineral::get_molar_volume)
    / map_endmembers_to_array(&Mineral::get_isothermal_bulk_modulus_reuss)
    * molar_fractions
  ).sum();
  return get_molar_volume()
    * 1.0 / (V_over_KT + solution_model->compute_VoverKT_excess());
}

double Solution::compute_isentropic_bulk_modulus_reuss() const {
  if (get_temperature() < constants::precision::abs_tolerance) {
    return get_isothermal_bulk_modulus_reuss();
  } else {
    return get_isothermal_bulk_modulus_reuss()
      * get_molar_heat_capacity_p()
      / get_molar_heat_capacity_v();
  }
}

double Solution::compute_isothermal_compressibility_reuss() const {
  return 1.0 / get_isothermal_bulk_modulus_reuss();
}

double Solution::compute_isentropic_compressibility_reuss() const {
  return 1.0 / get_isentropic_bulk_modulus_reuss();
}

double Solution::compute_shear_modulus() const {
  Eigen::ArrayXd em_shear_moduli = map_endmembers_to_array(&Mineral::get_shear_modulus);
  return averaging::utils::reuss_fn(molar_fractions, em_shear_moduli);
}

double Solution::compute_p_wave_velocity() const {
  constexpr double FOUR_THIRDS = 4.0 / 3.0;
  return std::sqrt(
    (get_isentropic_bulk_modulus_reuss() + FOUR_THIRDS * get_shear_modulus())
    / get_density()
  );
}

double Solution::compute_bulk_sound_velocity() const {
  return std::sqrt(get_isentropic_bulk_modulus_reuss() / get_density());
}

double Solution::compute_shear_wave_velocity() const {
  return std::sqrt(get_shear_modulus() / get_density());
}

double Solution::compute_grueneisen_parameter() const {
  if (get_temperature() < constants::precision::abs_tolerance) {
    return std::nan("");
  } else {
    return get_thermal_expansivity()
      * get_isothermal_bulk_modulus_reuss()
      * get_molar_volume()
      / get_molar_heat_capacity_v();
  }
}

double Solution::compute_thermal_expansivity() const {
  double alphaV = (
    map_endmembers_to_array(&Mineral::get_molar_volume)
    * map_endmembers_to_array(&Mineral::get_thermal_expansivity)
    * molar_fractions
  ).sum();
  return (1.0 / get_molar_volume())
    * (alphaV + solution_model->compute_alphaV_excess());
}

double Solution::compute_molar_heat_capacity_v() const {
  return get_molar_heat_capacity_p()
    - get_molar_volume()
    * get_temperature()
    * get_thermal_expansivity()
    * get_thermal_expansivity()
    * get_isothermal_bulk_modulus_reuss();
}

double Solution::compute_molar_heat_capacity_p() const {
  Eigen::ArrayXd em_Cp = map_endmembers_to_array(&Mineral::get_molar_heat_capacity_p);
  return (em_Cp * molar_fractions).sum() + solution_model->compute_Cp_excess();
}

double Solution::compute_isentropic_thermal_gradient() const {
  return
    (get_molar_volume() * get_temperature() * get_thermal_expansivity())
    / get_molar_heat_capacity_p();
}

types::FormulaMap Solution::compute_formula() const {
  return utils::sum_formulae(get_endmember_formulae(), molar_fractions);
}

// Setup functions for CompositeMaterial properties

int Solution::compute_n_endmembers() const {
  return static_cast<int>(solution_model->n_endmembers);
}

void Solution::setup_endmember_names() const {
  set_endmember_names(map_endmembers_to_vector<std::string>(&Mineral::get_name));
}

void Solution::setup_endmember_formulae() const {
  set_endmember_formulae(map_endmembers_to_vector<types::FormulaMap>(&Mineral::get_formula));
}

} // namespace burnman
