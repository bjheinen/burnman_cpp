/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/core/assemblage.hpp"
#include "burnman/core/mineral.hpp"
#include "burnman/utils/utils.hpp"

void Assemblage::reset() {
  // Reset caches Material properties
  Material::reset();
  // Reset cached Assemblage properties
  volume_fractions.reset();
}

// Public setters for assemblage properties
void Assemblage::set_averaging_scheme(AveragingType scheme_type) {
  switch (scheme_type) {
    case AveragingType::Voigt:
      averaging_scheme = std::make_unique<Voigt>();
      break;
    case AveragingType::Reuss:
      averaging_scheme = std::make_unique<Reuss>();
      break;
    case AveragingType::VoigtReussHill:
      averaging_scheme = std::make_unique<VoigtReussHill>();
      break;
    case AveragingType::HashinShtrikmanLower:
      averaging_scheme = std::make_unique<HashinShtrikmanLower>();
      break;
    case AveragingType::HashinShtrikmanUpper:
      averaging_scheme = std::make_unique<HashinShtrikmanUpper>();
      break;
    case AveragingType::HashinShtrikman:
      averaging_scheme = std::make_unique<HashinShtrikman>();
      break;
    default:
      throw std::invalid_argument("Unknown Averaging Scheme!");
  }
}

void Assemblage::set_averaging_scheme(std::unique_ptr<Averaging> custom_scheme) {
  averaging_scheme = std::move(custom_scheme);
}



// Public setter overrides of Material
void Assemblage::set_method(std::shared_ptr<EquationOfState> new_method) {
  for (auto& ph : phases) {
    ph->set_method(new_method);
  }
  // Clear properties cache
  reset();
}

void Assemblage::set_method(EOSType new_method) {
  for (auto& ph : phases) {
    ph->set_method(new_method);
  }
  // Clear properties cache
  reset();
}

void Assemblage::set_state(
  double new_pressure,
  double new_temperature
) {
  reset(); // TODO: check if reset needed here???
  // Set P,T using Material
  Material::set_state(new_pressure, new_temperature);
  // Set state of each endmember
  for (auto& ph : phases) {
    ph->set_state(new_pressure, new_temperature);
  }
}

// Public getter functions with caching
Eigen::ArrayXd Assemblage::get_volume_fractions() const {
  if (!volume_fractions.has_value()) {
    volume_fractions = compute_volume_fractions();
  }
  return *volume_fractions;
}

std::vector<int> Assemblage::get_endmembers_per_phase() const {
  if (!endmembers_per_phase.has_value()) {
    setup_endmember_properties();
  }
  return *endmembers_per_phase;
}

// Material property overrides
double Assemblage::compute_molar_internal_energy() const {
 return (map_phases_to_array(&Material::get_molar_internal_energy)
   * molar_fractions).sum();
}

double Assemblage::compute_molar_gibbs() const {
  return (map_phases_to_array(&Material::get_molar_gibbs)
    * molar_fractions).sum();
}

double Assemblage::compute_molar_helmholtz() const {
  return (map_phases_to_array(&Material::get_molar_helmholtz)
    * molar_fractions).sum();
}

double Assemblage::compute_molar_mass() const {
  return (map_phases_to_array(&Material::get_molar_mass)
    * molar_fractions).sum();
}

double Assemblage::compute_molar_volume() const {
  return (map_phases_to_array(&Material::get_molar_volume)
    * molar_fractions).sum();
}

double Assemblage::compute_density() const {
  Eigen::ArrayXd densities = map_phases_to_array(&Material::get_density);
  Eigen::ArrayXd volumes = get_volume_fractions();
  return averaging_scheme->average_density(volumes, densities);
}

double Assemblage::compute_molar_entropy() const {
  return (map_phases_to_array(&Material::get_molar_entropy)
    * molar_fractions).sum();
}

double Assemblage::compute_molar_enthalpy() const {
  return (map_phases_to_array(&Material::get_molar_enthalpy)
    * molar_fractions).sum();
}

double Assemblage::compute_isothermal_bulk_modulus_reuss() const {
  Eigen::ArrayXd V_frac = get_volume_fractions();
  Eigen::ArrayXd K_ph = map_phases_to_array(&Material::get_isothermal_bulk_modulus_reuss);
  Eigen::ArrayXd G_ph = map_phases_to_array(&Material::get_shear_modulus);
  return averaging_scheme->average_bulk_moduli(V_frac, K_ph, G_ph);
}

double Assemblage::compute_isentropic_bulk_modulus_reuss() const {
  Eigen::ArrayXd V_frac = get_volume_fractions();
  Eigen::ArrayXd K_ph = map_phases_to_array(&Material::get_isentropic_bulk_modulus_reuss);
  Eigen::ArrayXd G_ph = map_phases_to_array(&Material::get_shear_modulus);
  return averaging_scheme->average_bulk_moduli(V_frac, K_ph, G_ph);
}

double Assemblage::compute_isothermal_compressibility_reuss() const {
  return 1.0 / get_isothermal_bulk_modulus_reuss();
}

double Assemblage::compute_isentropic_compressibility_reuss() const {
  return 1.0 / get_isentropic_bulk_modulus_reuss();
}

double Assemblage::compute_shear_modulus() const {
  Eigen::ArrayXd V_frac = get_volume_fractions();
  Eigen::ArrayXd K_ph = map_phases_to_array(&Material::get_isentropic_bulk_modulus_reuss);
  Eigen::ArrayXd G_ph = map_phases_to_array(&Material::get_shear_modulus);
  return averaging_scheme->average_shear_moduli(V_frac, K_ph, G_ph);
}

// TODO: think about redunancy/replication with `Solution::compute_p_wave...'
double Assemblage::compute_p_wave_velocity() const {
  constexpr double FOUR_THIRDS = 4.0 / 3.0;
  return std::sqrt(
    (get_isentropic_bulk_modulus_reuss() + FOUR_THIRDS * get_shear_modulus())
    / get_density()
  );
}

double Assemblage::compute_bulk_sound_velocity() const {
  return std::sqrt(
    get_isentropic_bulk_modulus_reuss() / get_density()
  );
}

double Assemblage::compute_shear_wave_velocity() const {
  return std::sqrt(
    get_shear_modulus() / get_density()
  );
}

double Assemblage::compute_grueneisen_parameter() const {
  return get_thermal_expansivity()
    * get_isothermal_bulk_modulus_reuss()
    * get_molar_volume()
    / get_molar_heat_capacity_v();
}

double Assemblage::compute_thermal_expansivity() const {
  Eigen::ArrayXd volumes = get_volume_fractions();
  Eigen::ArrayXd alphas = map_phases_to_array(&Material::get_thermal_expansivity);
  return averaging_scheme->average_thermal_expansivity(volumes, alphas);
}

double Assemblage::compute_molar_heat_capacity_v() const {
  Eigen::ArrayXd c_v = map_phases_to_array(&Material::get_molar_heat_capacity_v);
  return averaging_scheme->average_heat_capacity_v(molar_fractions, c_v);
}

double Assemblage::compute_molar_heat_capacity_p() const {
  Eigen::ArrayXd c_p = map_phases_to_array(&Material::get_molar_heat_capacity_p);
  return averaging_scheme->average_heat_capacity_p(molar_fractions, c_p);
}

FormulaMap Assemblage::compute_formula() const {
  std::vector<FormulaMap> phase_formulae =
    map_phases_to_vector<FormulaMap>(&Material::get_formula);
  return utils::sum_formulae(phase_formulae, molar_fractions);
}

// Compute / setup functions for assemblage properties

Eigen::ArrayXd Assemblage::compute_volume_fractions() const {
  return map_phases_to_array(&Material::get_molar_volume) * molar_fractions;
}

int Assemblage::compute_n_endmembers() const {
  return static_cast<int>(get_endmember_names().size());
}

void Assemblage::setup_endmember_names() const {
  setup_endmember_properties();
}

void Assemblage::setup_endmember_formulae() const {
  setup_endmember_properties();
}

void Assemblage::setup_endmember_properties() const {
  std::vector<FormulaMap> temp_formulae;
  std::vector<std::string> temp_names;
  std::vector<int> temp_empp;
  // Loop through phases
  for (const auto& ph : phases) {
    // If phase is a Solution
    if (auto sol = std::dynamic_pointer_cast<CompositeMaterial>(ph)) {
      // Get properties
      const auto& names = sol->get_endmember_names();
      const auto& formulae = sol->get_endmember_formulae();
      // Batch add endmember formulas
      temp_formulae.insert(temp_formulae.end(),
                                formulae.begin(), formulae.end());
      // Add endmember names
      temp_names.reserve(temp_names.size() + names.size());
      const std::string& phase_name  = sol->get_name();
      for (const auto& name_i : names) {
        temp_names.push_back(name_i + " in " + phase_name);
      }
      // Add number of endmembers to list
      temp_empp.push_back(sol->get_n_endmembers());
    // If phase is a Mineral
    } else if (auto min = std::dynamic_pointer_cast<Mineral>(ph)) {
      // Store mineral properties
      temp_formulae.push_back(min->get_formula());
      temp_names.push_back(min->get_name());
      temp_empp.push_back(1);
    } else {
      throw std::runtime_error(
        "Unsupported Material type; expected Mineral or Solution.");
    }
  }
  // Set properties
  set_endmember_names(temp_names);
  set_endmember_formulae(temp_formulae);
  set_endmembers_per_phase(temp_empp);
}

void Assemblage::set_endmembers_per_phase(std::vector<int> v) const {
  endmembers_per_phase = std::move(v);
}

/*
stoichiometric_array
reaction_basis_as_strings (skip for now)
reduced_stoichiometric_array
*/
