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
#include <cmath>
#include <stdexcept>
#include <string>
#include "burnman/utils/chemistry_utils.hpp"
#include "burnman/core/mineral.hpp"

namespace burnman {

void Assemblage::reset_cache() {
  // Reset cached Material & CompositeMaterial properties
  CompositeMaterial::reset_cache();
  // Reset cached Assemblage properties
  volume_fractions.reset();
  reaction_affinities.reset();
  equilibrium_tolerance = 1.0e-3;
}

void Assemblage::clear_computed_properties() {
  // Reset CompositeMaterial properties
  // CompositeMaterial::clear_computed_properties calls
  // Material::clear_computed_properties, which calls
  // reset_cache - virtual dispatch should result in
  // Assemblage::reset_cache being called.
  CompositeMaterial::clear_computed_properties();
  this->endmembers_per_phase.reset();
  // Clear properties stored by phases
  for (auto& ph : phases) {
    ph->clear_computed_properties();
  }
}

// Public setters for assemblage properties

// Single phase by ptr
void Assemblage::add_phase(const std::shared_ptr<Material>& phase_ptr) {
  this->phases.push_back(phase_ptr);
}

// Multiple phases by pointer
void Assemblage::add_phases(
  const std::vector<std::shared_ptr<Material>>& phase_ptr_list
) {
  this->phases.insert(
    this->phases.end(),
    phase_ptr_list.begin(),
    phase_ptr_list.end()
  );
}

void Assemblage::add_phases(
  std::initializer_list<std::shared_ptr<Material>> phase_ptr_list
) {
  this->phases.insert(
    this->phases.end(),
    phase_ptr_list.begin(),
    phase_ptr_list.end()
  );
}

void Assemblage::set_averaging_scheme(types::AveragingType scheme_type) {
  switch (scheme_type) {
    case types::AveragingType::Voigt:
      this->averaging_scheme = std::make_shared<averaging::Voigt>();
      break;
    case types::AveragingType::Reuss:
      this->averaging_scheme = std::make_shared<averaging::Reuss>();
      break;
    case types::AveragingType::VRH:
      this->averaging_scheme = std::make_shared<averaging::VoigtReussHill>();
      break;
    case types::AveragingType::HashinShtrikmanLower:
      this->averaging_scheme = std::make_shared<averaging::HashinShtrikmanLower>();
      break;
    case types::AveragingType::HashinShtrikmanUpper:
      this->averaging_scheme = std::make_shared<averaging::HashinShtrikmanUpper>();
      break;
    case types::AveragingType::HashinShtrikman:
      this->averaging_scheme = std::make_shared<averaging::HashinShtrikman>();
      break;
    default:
      throw std::invalid_argument("Unknown Averaging Scheme!");
  }
}

void Assemblage::set_averaging_scheme(std::shared_ptr<averaging::AveragingScheme> custom_scheme) {
  this->averaging_scheme = std::move(custom_scheme);
}

// Eigen::ArrayXd version (core implementation)
void Assemblage::set_fractions(
  const Eigen::ArrayXd& fractions,
  const types::FractionType fraction_type
) {
  // Asserts for valid fraction array
  if (fractions.size() != static_cast<Eigen::Index>(this->phases.size())) {
    throw std::invalid_argument(
      "Fractions size doesn't match number of phases!"
    );
  }
  if ((fractions < -1e-12).any()) {
    throw std::invalid_argument(
      "Fractions contain negative values! (Less than -1e-12)."
    );
  }
  double total = fractions.sum();
  reset_cache();
  Eigen::ArrayXd norm_fractions = fractions;
  if (std::abs(total - 1.0) > 1e-12) {
    // Todo: Warnings"
    //std::cerr << "Warning: list of fractions does not add up to one but "
    //          << total << ". Normalizing.\n";
    norm_fractions /= total;
  }
  switch (fraction_type) {
    case types::FractionType::Molar :
      this->molar_fractions = norm_fractions;
      break;
    case types::FractionType::Mass :
      this->molar_fractions = convert_mass_to_molar_fractions(norm_fractions);
      break;
    default :
      throw std::invalid_argument(
        "Fraction type not recognised. Please use 'molar' or 'mass'."
      );
  }
  // Clip to zero as in Py version
  this->molar_fractions = molar_fractions.cwiseMax(0.0);
}

void Assemblage::set_fractions(
  std::initializer_list<double> fractions,
  const types::FractionType fraction_type
) {
  this->set_fractions(
    Eigen::Map<const Eigen::ArrayXd>(fractions.begin(), fractions.size()),
    fraction_type
  );
}

void Assemblage::set_n_moles(double new_n_moles) {
  this->n_moles = new_n_moles;
}

void Assemblage::set_equilibrium_tolerance(double new_equilibrium_tolerance) {
  this->equilibrium_tolerance = new_equilibrium_tolerance;
}

// Public setter overrides of Material
void Assemblage::set_method(std::shared_ptr<EquationOfState> new_method) {
  for (auto& ph : phases) {
    ph->set_method(new_method);
  }
  // Clear properties cache
  reset_cache();
}

void Assemblage::set_method(types::EOSType new_method) {
  for (auto& ph : phases) {
    ph->set_method(new_method);
  }
  // Clear properties cache
  reset_cache();
}

void Assemblage::set_state(
  double new_pressure,
  double new_temperature
) {
  reset_cache(); // TODO: check if reset needed here???
  // Set P,T using Material
  Material::set_state(new_pressure, new_temperature);
  // Set state of each endmember
  for (auto& ph : phases) {
    ph->set_state(new_pressure, new_temperature);
  }
}

// Public convenience getters
Eigen::Index Assemblage::get_n_phases() const {
  return static_cast<Eigen::Index>(this->phases.size());
}

std::shared_ptr<Material> Assemblage::get_phase(std::size_t index) const {
  return this->phases.at(index);
}

Eigen::ArrayXd Assemblage::get_molar_fractions() const {
  return this->molar_fractions;
}

double Assemblage::get_n_moles() const {
  if (!this->n_moles.has_value()) {
    throw std::runtime_error("n_moles not set!");
  }
  return *n_moles;
}

double Assemblage::get_equilibrium_tolerance() const {
  return this->equilibrium_tolerance;
}

// Public getter functions with caching
Eigen::ArrayXd Assemblage::get_volume_fractions() const {
  if (!volume_fractions.has_value()) {
    volume_fractions = compute_volume_fractions();
  }
  return *volume_fractions;
}

Eigen::VectorXd Assemblage::get_reaction_affinities() const {
  if (!reaction_affinities.has_value()) {
    reaction_affinities = compute_reaction_affinities();
  }
  return *reaction_affinities;
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

types::FormulaMap Assemblage::compute_formula() const {
  std::vector<types::FormulaMap> phase_formulae =
    map_phases_to_vector<types::FormulaMap>(&Material::get_formula);
  return utils::sum_formulae(phase_formulae, molar_fractions);
}

// Compute / setup functions for assemblage properties

Eigen::ArrayXd Assemblage::compute_volume_fractions() const {
  return map_phases_to_array(&Material::get_molar_volume) * molar_fractions;
}

Eigen::ArrayXd Assemblage::compute_partial_gibbs() const {
  Eigen::ArrayXd embr_partial_gibbs(get_n_endmembers());
  // Loop over phases
  Eigen::Index j = 0;
  for (const auto& ph : this->phases) {
    Eigen::Index n = 1;
    // If phase is a Solution or Assemblage
    if (auto sol = std::dynamic_pointer_cast<CompositeMaterial>(ph)) {
      n = sol->get_n_endmembers();
      embr_partial_gibbs.segment(j, n) = sol->get_partial_gibbs();
    } else {
      embr_partial_gibbs(j) = ph->get_molar_gibbs();
    }
    j += n;
  }
  return embr_partial_gibbs;
}

Eigen::VectorXd Assemblage::compute_reaction_affinities() const {
  return get_reaction_basis() * get_partial_gibbs().matrix();
}

Eigen::Index Assemblage::compute_n_endmembers() const {
  return static_cast<Eigen::Index>(get_endmember_names().size());
}

void Assemblage::setup_endmember_names() const {
  setup_endmember_properties();
}

void Assemblage::setup_endmember_formulae() const {
  setup_endmember_properties();
}

void Assemblage::setup_endmember_properties() const {
  std::vector<types::FormulaMap> temp_formulae;
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
      temp_empp.push_back(static_cast<int>(sol->get_n_endmembers()));
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

Eigen::ArrayXd Assemblage::convert_mass_to_molar_fractions(
  const Eigen::ArrayXd& mass_fractions
) const {
  Eigen::ArrayXd moles = mass_fractions / map_phases_to_array(&Material::get_molar_mass);
  return moles / moles.sum();
}

} // namespace burnman
