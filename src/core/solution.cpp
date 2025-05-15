/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <cmath>
#include "burnman/core/solution.hpp"
#include "burnman/utils/constants.hpp"
#include "burnman/utils/chemistry_utils.hpp"
#include "burnman/utils/matrix_utils.hpp"
#include "burnman/core/averaging_schemes.hpp"

void Solution::reset() {
  // Reset caches Material properties
  Material::reset();
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
  partial_gibbs.reset();
  partial_volumes.reset();
  partial_entropies.reset();
  gibbs_hessian.reset();
  entropy_hessian.reset();
  volume_hessian.reset();
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

Eigen::ArrayXd Solution::get_partial_gibbs() const {
  if (!partial_gibbs.has_value()) {
    partial_gibbs = compute_partial_gibbs();
  }
  return *partial_gibbs;
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
  return averaging::reuss_fn(molar_fractions, em_shear_moduli);
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

// Setup functions for solution properties stored on class initialisation

void Solution::setup_endmember_names() {
  endmember_names = map_endmembers_to_vector<std::string>(&Mineral::get_name);
}

void Solution::setup_endmember_formulae() {
  endmember_formulae = map_endmembers_to_vector<FormulaMap>(&Mineral::get_formula);
}

void Solution::setup_elements() {
  // Grab set of all elements in formulae
  std::unordered_set<std::string> all_elements;
  for (const auto& embr : endmember_formulae) {
    for (const auto& [element, n] : embr) {
      all_elements.insert(element);
    }
  }
  // Sort elements into IUPAC order
  elements = utils::sort_element_list_to_IUPAC_order(all_elements);
}

void Solution::setup_stoichiometric_matrix() {
  int n_elements = static_cast<int>(elements.size());
  int n_endmembers = static_cast<int>(solution_model->n_endmembers);
  stoichiometric_matrix.resize(n_endmembers, n_elements);
  stoichiometric_matrix.setZero();
  for (int i = 0; i < n_endmembers; ++i) {
    const auto& formula = endmember_formulae[i];
    for (int j = 0; j < n_elements; ++j) {
      const std::string& element = elements[j];
      auto it = formula.find(element);
      if (it != formula.end()) {
        stoichiometric_matrix(i, j) = it->second;
      }
    }
  }
}

void Solution::setup_independent_element_indices() {
  independent_element_indices.clear();
  independent_element_indices = utils::get_independent_row_indices(stoichiometric_matrix);
}

void Solution::setup_dependent_element_indices() {
  dependent_element_indices.clear();
  for (int i = 0; i < elements.size(); ++i) {
    if (std::find(
      independent_element_indices.begin(),
      independent_element_indices.end(), i)
        == independent_element_indices.end()
    ) {
      dependent_element_indices.push_back(i);
    }
  }
}

void Solution::setup_reaction_basis() {
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(stoichiometric_matrix.transpose());
  Eigen::MatrixXd nullspace = lu_decomp.kernel();
  // Maybe consider threshold?
  if (nullspace.cols() == 0) {
    return Eigen::MatrixXd(0, solution_model->n_endmembers);
  }
  return nullspace.transpose();
}

void Solution::setup_n_reactions() {
  n_reactions = reaction_basis.rows();
}

void Solution::setup_compositional_basis() {

}

void Solution::setup_compositional_null_basis() {

}

void Solution::setup_solution_properties() {

}
