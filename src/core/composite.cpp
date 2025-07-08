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
#include "burnman/utils/utils.hpp"
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

int Composite::get_n_elements() const {
  if (!n_elements.has_value()) {
    setup_endmember_properties();
  }
  return *n_elements;
}

int Composite::get_n_endmembers() const {
  if (!n_endmembers.has_value()) {
    setup_endmember_properties();
  }
  return *n_endmembers;
}

std::vector<int> Composite::get_endmembers_per_phase() const {
  if (!endmembers_per_phase.has_value()) {
    setup_endmember_properties();
  }
  return *endmembers_per_phase;
}

std::vector<std::string> Composite::get_elements() const {
  if (!elements.has_value()) {
    setup_endmember_properties();
  }
  return *elements;
}

std::vector<std::string> Composite::get_endmember_names() const {
  if (!endmember_names.has_value()) {
    setup_endmember_properties();
  }
  return *endmember_names;
}

std::vector<FormulaMap> Composite::get_endmember_formulae() const {
  if (!endmember_formulae.has_value()) {
    setup_endmember_properties();
  }
  return *endmember_formulae;
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

// Compute / setup functions for composite properties

Eigen::ArrayXd Composite::compute_volume_fractions() const {
  return map_phases_to_array(&Material::get_molar_volume) * molar_fractions;
}

void Composite::setup_endmember_properties() const {

  // This should only be called once,
  // but check for value and reset just in case

  if (endmember_formulae.has_value()) {
    endmember_formulae->clear();
  } else {
    endmember_formulae = std::vector<FormulaMap>{};
  }

  if (endmember_names.has_value()) {
    endmember_names->clear();
  } else {
    endmember_names = std::vector<std::string>{};
  }

  if (endmembers_per_phase.has_value()) {
    endmembers_per_phase->clear();
  } else {
    endmembers_per_phase = std::vector<int>{};
  }

  if (elements.has_value()) {
    elements->clear();
  } else {
    elements = std::vector<std::string>{};
  }

  // Loop through phases
  for (const auto& ph : phases) {
    // If phase is a Solution
    if (auto sol = std::dynamic_pointer_cast<Solution>(ph)) {
      // Get properties
      const auto& names = sol->get_endmember_names();
      const auto& formulae = sol->get_endmember_formulae();
      // Batch add endmember formulas
      endmember_formulae->insert(endmember_formulae->end(),
                                 formulae.begin(), formulae.end());
      // Add endmember names
      endmember_names->reserve(endmember_names->size() + names.size());
      const std::string& phase_name  = sol->get_name();
      for (const auto& name : names) {
        endmember_names->push_back(name + " in " + phase_name);
      }
      // Add number of endmembers to list
      endmembers_per_phase->push_back(
        static_cast<int>(sol->get_n_endmembers()));
    // If phase is a Mineral
    } else if (auto min = std::dynamic_pointer_cast<Mineral>(ph)) {
      // Store mineral properties
      endmember_formulae->push_back(min->get_formula());
      endmember_names->push_back(min->get_name());
      endmembers_per_phase->push_back(1);
    } else {
      throw std::runtime_error(
        "Unsupported Material type; expected Mineral or Solution.");
    }
  }
  // Grab set of all elements in formulae
  std::unordered_set<std::string> all_elements;
  for (const auto& embr : *endmember_formulae) {
    for (const auto& [element, n] : embr) {
      all_elements.insert(element);
    }
  }
  // Sort elements into IUPAC order and set property
  elements = utils::sort_element_list_to_IUPAC_order(all_elements);
  // Set derived properties
  n_elements = elements->size();
  n_endmembers = endmember_names->size();
}

/*
stoichiometric_matrix
stoichiometric_array
reaction_basis
reaction_basis_as_strings (skip for now)
n_reactions
independent_element_indices
dependent_element_indices
reduced_stoichiometric_array
compositional_null_basis
*/
