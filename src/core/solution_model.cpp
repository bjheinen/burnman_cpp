/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <algorithm> // For std::count
#include <regex>
#include <stdexcept>

#include "burnman/core/solution_model.hpp"
#include "burnman/utils/constants.hpp"
#include "burnman/utils/string_utils.hpp"
#include "burnman/utils/math_utils.hpp"
#include "burnman/utils/matrix_utils.hpp"

SolutionModel::SolutionModel(const PairedEndmemberList& endmember_list) {
  n_endmembers = endmember_list.size();
  endmembers.reserve(n_endmembers);
  formulas.reserve(n_endmembers);
  // Unpack and store endmembers/formulas separately
  for (const auto& [mineral, formula] : endmember_list) {
    endmembers.push_back(mineral);
    formulas.push_back(formula);
  }
  // Process chemistry sets up site multiplicities/occupancies etc.
  process_solution_chemistry();
}

void SolutionModel::process_solution_chemistry() {

  // Reset n_occupancies count to zero
  n_occupancies = 0;

  // Set class n_sites
  n_sites = std::count(formulas[0].begin(), formulas[0].end(), '[');

  // Check that number of sites is constant
  for (const std::string& f : formulas) {
    if (std::count(f.begin(), f.end(), '[') != n_sites) {
      throw std::runtime_error("All formulae must have the same number of distinct sites.");
    }
  }

  // Make intermediate occ/mult data containers for parsing
  Eigen::ArrayXXd formula_multiplicities;
  std::vector<std::vector<std::vector<double>>> list_occupancies;

  // Store multiplicities
  formula_multiplicities.resize(n_endmembers, n_sites);
  // Resize other data containers
  solution_formulae.resize(n_endmembers);
  sites.resize(n_sites);
  // List occupancies is triple nested vector
  // std::vector<std::vector<std::vector<double>>>
  // Outer --> n_endmembers
  // Inner --> n_site
  // Innermost --> n_species (not known until parsed)
  list_occupancies.resize(n_endmembers);
  for (int i = 0; i < n_endmembers; ++i) {
    list_occupancies[i].resize(n_sites);
    for (int j = 0; j < n_sites; ++j) {
      list_occupancies[i][j].resize(0);
    }
  }

  // TODO: move these to more senible place, maybe factor out into string_utils
  // Regex strings for splits
  std::regex site_split_regex(R"(\[)");
  std::regex occ_split_regex(R"(\])");
  std::regex species_split_regex("[A-Z][^A-Z]*");
  std::regex species_frac_split_regex("([0-9][^A-Z]*)");

  // Loop over endmembers
  for (int i_mbr = 0; i_mbr < n_endmembers; ++i_mbr) {
    // Split formula into sites - 'Mg]3', 'Al]2', 'etc.'
    std::sregex_token_iterator it_sites(
      formulas[i_mbr].begin(), formulas[i_mbr].end(),
      site_split_regex, -1);
    // Discard string before first [ by ++it first
    std::vector<std::string> site_formulas(++it_sites, {});
    // Loop over sites in formula
    for (int i_site = 0; i_site < n_sites; ++i_site) {
      // Split on ] to get site occupancy and multiplicity
      std::sregex_token_iterator it_occ(
        site_formulas[i_site].begin(), site_formulas[i_site].end(),
        occ_split_regex, -1);
      // To check if match exists before dereferencing,
      // do if (it != std::sregex_token_iterator())
      std::string site_occ = *it_occ++;
      std::string site_mult_str = (it_occ != std::sregex_token_iterator()) ? it_occ->str() : "";

      // site multiplicity may have rest of formula, so split
      site_mult_str = utils::extract_numeric_prefix(site_mult_str);
      double site_mult = site_mult_str.empty() ? 1.0 : utils::stod(site_mult_str);

      // Store multiplicity
      formula_multiplicities(i_mbr, i_site) = site_mult;

      // Split occupancy into species
      std::vector<std::string> species;
      std::sregex_iterator it_species(site_occ.begin(), site_occ.end(), species_split_regex);
      for (std::sregex_iterator i = it_species; i != std::sregex_iterator(); ++i) {
        species.push_back(i->str());
      }

      // Loop over species on site
      for (const std::string& sp : species) {

        // "Mg" --> "Mg"; "Mg1/2" --> ["Mg", "1/2"]
        std::sregex_token_iterator it_frac(
          sp.begin(), sp.end(),
          species_frac_split_regex, {-1, 1}
        );
        std::vector<std::string> species_split(it_frac, {});
        // Get species name and proportion (name = Mg, prop = 0.5 etc.)
        std::string species_name = species_split[0];
        double proportion_species_on_site;
        if (species_split.size() == 1) {
          proportion_species_on_site = 1.0;
        } else {
          proportion_species_on_site = utils::stod(species_split[1]);
        }

        // TODO: When species spread across sites?
        solution_formulae[i_mbr][species_name] += site_mult * proportion_species_on_site;

        // Add to sites[i_site] if not present
        auto& site_species = sites[i_site];
        auto species_pos = std::find(site_species.begin(), site_species.end(), species_name);
        int i_el;
        if (species_pos == site_species.end()) {
          // Use current size as index of next entry
          i_el = site_species.size();
          site_species.push_back(species_name);
          ++n_occupancies;
          // Append 0 to list_occupancies for already parsed endmembers
          // when found new species
          for (int k = 0; k <= i_mbr; ++k) {
            // Check and resize outer if needed
            auto& occupancies = list_occupancies[k];
            if (occupancies.size() <= i_site) {
              occupancies.resize(i_site + 1);
            }
            // Check and resize inner if needed
            auto& site_occupancy = occupancies[i_site];
            if (site_occupancy.size() <= i_el) {
              site_occupancy.resize(i_el + 1, 0.0);
            }
          }
        } else {
          i_el = std::distance(site_species.begin(), species_pos);
          // Check and resize inner if needed
          auto& site_occupancy = list_occupancies[i_mbr][i_site];
          if (site_occupancy.size() <= i_el) {
            site_occupancy.resize(i_el + 1, 0.0);
          }
        }
        // Store species occupancy
        list_occupancies[i_mbr][i_site][i_el] = proportion_species_on_site;
      } // Species loop
      // Could add here looping over species not on site to add to solution_formulae map
      // Solution formulae currently unused so ignoring for now
    } // Site loop
  } // Endmember loop

  // Resize occupancy/multiplicity arrays
  endmember_occupancies.resize(n_endmembers, n_occupancies);
  site_multiplicities.resize(n_endmembers, n_occupancies);

  // Loop over endmembers again
  for (int i_mbr = 0; i_mbr < n_endmembers; ++i_mbr) {
    int n_species = 0;
    for (int i_site = 0; i_site < n_sites; ++i_site) {
      for (size_t i_el = 0; i_el < list_occupancies[i_mbr][i_site].size(); ++i_el) {
        endmember_occupancies(i_mbr, n_species) = list_occupancies[i_mbr][i_site][i_el];
        site_multiplicities(i_mbr, n_species) = formula_multiplicities(i_mbr, i_site);
        ++n_species;
      }
    }
  } // Endmember loop

  // Element-wise multiply for n_occ
  endmember_n_occupancies = endmember_occupancies * site_multiplicities;

  // Get site names
  site_names.clear();
  for (int i_site = 0; i_site < n_sites; ++i_site) {
    // Grab uppercase letter from index (works to 26!)
    char site_id = 'A' + i_site;
    for (const std::string& sp : sites[i_site]) {
      site_names.push_back(sp + "_" + site_id);
    }
  }

  // Do parsed chemical formula etc. (general and empty)
  // Replace [sites] with [] on one formula for empty formula
  empty_formula = std::regex_replace(formulas[0], std::regex(R"(\[.*?\])"), "[]");
  // Split empty formula on [, then re-join with site species for general
  std::regex split_brackets(R"(\[)");
  std::sregex_token_iterator it_empty_split(
    empty_formula.begin(),
    empty_formula.end(),
    split_brackets, -1);
  std::sregex_token_iterator end;
  std::vector<std::string> split_empty(it_empty_split, end);
  // Cut first element here if needed to stop overflow from pre-increment
  if (!split_empty.empty() && split_empty[0].empty())
    split_empty.erase(split_empty.begin());
  general_formula = split_empty[0];
  for (int i = 0; i < n_sites; ++i) {
    general_formula += "[" + utils::join(sites[i], ",") + split_empty[i + 1];
  }

}

double SolutionModel::compute_excess_gibbs_free_energy(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return (
    molar_fractions
    * compute_excess_partial_gibbs_free_energies(pressure, temperature, molar_fractions)
  ).sum();
}

double SolutionModel::compute_excess_volume(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return (
    molar_fractions
    * compute_excess_partial_volumes(pressure, temperature, molar_fractions)
  ).sum();
}

double SolutionModel::compute_excess_entropy(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return (
    molar_fractions
    * compute_excess_partial_entropies(pressure, temperature, molar_fractions)
  ).sum();
}

double SolutionModel::compute_excess_enthalpy(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return compute_excess_gibbs_free_energy(
      pressure, temperature, molar_fractions)
    + temperature * compute_excess_entropy(
      pressure, temperature, molar_fractions);
}

double SolutionModel::compute_Cp_excess() const {
  return 0.0;
}

double SolutionModel::compute_alphaV_excess() const {
  return 0.0;
}

double SolutionModel::compute_VoverKT_excess() const {
  return 0.0;
}

// Constructor for IdealSolution
IdealSolution::IdealSolution(const SolutionModel::PairedEndmemberList& endmember_list)
  : SolutionModel(endmember_list) {
  // Calculate configurational entropies also
  endmember_configurational_entropies = compute_endmember_configurational_entropies();
}

// Public function overrides for IdealSolution
Eigen::ArrayXd IdealSolution::compute_excess_partial_gibbs_free_energies(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return compute_ideal_excess_partial_gibbs(temperature, molar_fractions);
}

Eigen::ArrayXd IdealSolution::compute_excess_partial_entropies(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return compute_ideal_excess_partial_entropies(molar_fractions);
}

Eigen::ArrayXd IdealSolution::compute_excess_partial_volumes(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return Eigen::ArrayXd::Zero(n_endmembers);
}

Eigen::MatrixXd IdealSolution::compute_gibbs_hessian(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return -temperature * compute_ideal_entropy_hessian(molar_fractions);
}

Eigen::MatrixXd IdealSolution::compute_entropy_hessian(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return compute_ideal_entropy_hessian(molar_fractions);
}

Eigen::MatrixXd IdealSolution::compute_volume_hessian(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return Eigen::MatrixXd::Zero(n_endmembers, n_endmembers);
}

Eigen::ArrayXd IdealSolution::compute_activities(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return compute_ideal_activities(molar_fractions);
}

Eigen::ArrayXd IdealSolution::compute_activity_coefficients(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return Eigen::ArrayXd::Ones(n_endmembers);
}

// Private functions for IdealSolution
Eigen::ArrayXd IdealSolution::compute_endmember_configurational_entropies() const {
  return -constants::physics::gas_constant
  * (endmember_n_occupancies
    * (utils::logish(endmember_n_occupancies)
      - utils::logish(site_multiplicities))
    ).rowwise().sum();
}

Eigen::ArrayXd IdealSolution::compute_ideal_excess_partial_gibbs(
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return -temperature * compute_ideal_excess_partial_entropies(molar_fractions);
}

Eigen::ArrayXd IdealSolution::compute_ideal_excess_partial_entropies(
  const Eigen::ArrayXd& molar_fractions
) const {
  return -constants::physics::gas_constant
    * compute_log_ideal_activities(molar_fractions);
}

Eigen::ArrayXd IdealSolution::compute_ideal_activities(
  const Eigen::ArrayXd& molar_fractions
) const {
  // Dot product
  Eigen::ArrayXd reduced_n_occupancies = (endmember_n_occupancies.colwise() * molar_fractions).colwise().sum();
  Eigen::ArrayXd reduced_multiplicities = (site_multiplicities.colwise() * molar_fractions).colwise().sum();
  Eigen::ArrayXd reduced_occupancies = reduced_n_occupancies * utils::inverseish(reduced_multiplicities);
  // endmember_n_occupancies is class member
  double a = reduced_occupancies.pow(endmember_n_occupancies).prod();
  Eigen::ArrayXd norm_constants = (endmember_configurational_entropies / constants::physics::gas_constant).exp();
  return norm_constants * a;
}

Eigen::ArrayXd IdealSolution::compute_log_ideal_activities(
  const Eigen::ArrayXd& molar_fractions
) const {
  Eigen::ArrayXd reduced_n_occupancies = (endmember_n_occupancies.colwise() * molar_fractions).colwise().sum();
  Eigen::ArrayXd reduced_multiplicities = (site_multiplicities.colwise() * molar_fractions).colwise().sum();
  Eigen::ArrayXd lna = (
    endmember_n_occupancies.rowwise()
    * (utils::logish(reduced_n_occupancies)
      - utils::logish(reduced_multiplicities)
      ).transpose()
  ).rowwise().sum();
  Eigen::ArrayXd norm_constants = endmember_configurational_entropies / constants::physics::gas_constant;
  return lna + norm_constants;
}

Eigen::MatrixXd IdealSolution::compute_log_ideal_activity_derivatives(
  const Eigen::ArrayXd& molar_fractions
) const {
  Eigen::ArrayXd reduced_n_occupancies = (endmember_n_occupancies.colwise() * molar_fractions).colwise().sum();
  Eigen::ArrayXd reduced_multiplicities = (site_multiplicities.colwise() * molar_fractions).colwise().sum();
  Eigen::MatrixXd dlnadp =
    ((endmember_n_occupancies.rowwise() * utils::inverseish(reduced_n_occupancies).transpose()).matrix()
      * endmember_n_occupancies.matrix().transpose())
    - ((endmember_n_occupancies.rowwise() * utils::inverseish(reduced_multiplicities).transpose()).matrix()
      * site_multiplicities.matrix().transpose());
  return dlnadp;
}

Eigen::MatrixXd IdealSolution::compute_ideal_entropy_hessian(
  const Eigen::ArrayXd& molar_fractions
) const {
  return constants::physics::gas_constant
    * compute_log_ideal_activity_derivatives(molar_fractions);
}

// Constructor for AsymmetricRegularSolution
AsymmetricRegularSolution::AsymmetricRegularSolution(
  const SolutionModel::PairedEndmemberList& endmember_list,
  std::vector<double> alphas_vector,
  std::vector<std::vector<double>> energy_interaction,
  std::vector<std::vector<double>> volume_interaction,
  std::vector<std::vector<double>> entropy_interaction
) : IdealSolution(endmember_list) {
  // Map alphas to an Eigen::ArrayXd
  Eigen::ArrayXd alphas = Eigen::Map<Eigen::ArrayXd>(
    alphas_vector.data(), alphas_vector.size());

  W_e = utils::populate_interaction_matrix(
    utils::jagged2square(energy_interaction, n_endmembers),
    alphas,
    n_endmembers);

  if (!volume_interaction.empty()) {
    W_v = utils::populate_interaction_matrix(
      utils::jagged2square(volume_interaction, n_endmembers),
      alphas,
      n_endmembers);
  } else {
    W_v = Eigen::MatrixXd::Zero(n_endmembers, n_endmembers);
  }
  if (!entropy_interaction.empty()) {
    W_s = utils::populate_interaction_matrix(
      utils::jagged2square(entropy_interaction, n_endmembers),
      alphas,
      n_endmembers);
  } else {
    W_s = Eigen::MatrixXd::Zero(n_endmembers, n_endmembers);
  }
}

// Public function overrides for AsymmetricRegularSolution
Eigen::ArrayXd AsymmetricRegularSolution::compute_excess_partial_gibbs_free_energies(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return IdealSolution::compute_excess_partial_gibbs_free_energies(
      pressure, temperature, molar_fractions)
    + compute_non_ideal_excess_partial_gibbs(
      pressure, temperature, molar_fractions);
}

Eigen::ArrayXd AsymmetricRegularSolution::compute_excess_partial_entropies(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return IdealSolution::compute_excess_partial_entropies(
      pressure, temperature, molar_fractions)
    + compute_non_ideal_interactions(W_e, molar_fractions);
}

Eigen::ArrayXd AsymmetricRegularSolution::compute_excess_partial_volumes(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return IdealSolution::compute_excess_partial_volumes(
      pressure, temperature, molar_fractions)
    + compute_non_ideal_interactions(W_v, molar_fractions);
}

Eigen::MatrixXd AsymmetricRegularSolution::compute_gibbs_hessian(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  Eigen::MatrixXd interactions = W_e - temperature * W_s + pressure * W_v;
  return IdealSolution::compute_gibbs_hessian(
      pressure, temperature, molar_fractions)
    + compute_non_ideal_hessian(interactions, molar_fractions);
}

Eigen::MatrixXd AsymmetricRegularSolution::compute_entropy_hessian(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return IdealSolution::compute_entropy_hessian(
      pressure, temperature, molar_fractions)
    + compute_non_ideal_hessian(W_s, molar_fractions);
}

Eigen::MatrixXd AsymmetricRegularSolution::compute_volume_hessian(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return IdealSolution::compute_volume_hessian(
      pressure, temperature, molar_fractions)
    + compute_non_ideal_hessian(W_v, molar_fractions);
}

Eigen::ArrayXd AsymmetricRegularSolution::compute_activities(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return IdealSolution::compute_activities(
      pressure, temperature, molar_fractions)
    * compute_activity_coefficients(
      pressure, temperature, molar_fractions);
}

Eigen::ArrayXd AsymmetricRegularSolution::compute_activity_coefficients(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  if (temperature < constants::precision::abs_tolerance) {
    throw std::runtime_error("Activity coefficients undefined at 0 K.");
  }
  return (
    compute_non_ideal_excess_partial_gibbs(
      pressure, temperature, molar_fractions)
    / (constants::physics::gas_constant * temperature)
  ).exp();
}

// Private compute functions for AsymmetricRegularSolution

Eigen::ArrayXd AsymmetricRegularSolution::compute_phi(
  const Eigen::ArrayXd& molar_fractions
) const {
  Eigen::ArrayXd phi = alphas * molar_fractions;
  return phi / phi.sum();
}

Eigen::ArrayXd AsymmetricRegularSolution::compute_non_ideal_interactions(
  const Eigen::MatrixXd& W,
  const Eigen::ArrayXd& molar_fractions
) const {
  Eigen::ArrayXd phi = compute_phi(molar_fractions);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_endmembers, n_endmembers);
  Eigen::MatrixXd q = I.rowwise() - phi.matrix().transpose();
  Eigen::ArrayXd W_int = -alphas * (q * W).cwiseProduct(q).rowwise().sum().array();
  return W_int;
}

Eigen::ArrayXd AsymmetricRegularSolution::compute_non_ideal_excess_partial_gibbs(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  Eigen::ArrayXd E_int = compute_non_ideal_interactions(W_e, molar_fractions);
  Eigen::ArrayXd S_int = compute_non_ideal_interactions(W_s, molar_fractions);
  Eigen::ArrayXd V_int = compute_non_ideal_interactions(W_v, molar_fractions);
  return E_int - temperature * S_int + pressure * V_int;
}

Eigen::MatrixXd AsymmetricRegularSolution::compute_non_ideal_hessian(
  const Eigen::MatrixXd& interactions,
  const Eigen::ArrayXd& molar_fractions
) const {
  // Maybe factor out - reused in non_ideal_interactions
  Eigen::ArrayXd phi = compute_phi(molar_fractions);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_endmembers, n_endmembers);
  Eigen::MatrixXd q = I.rowwise() - phi.matrix().transpose();
  //
  double sum_pa = molar_fractions.matrix().dot(alphas.matrix());
  Eigen::MatrixXd alpha_outer_product = alphas.matrix() * (alphas/sum_pa).matrix().transpose();
  Eigen::MatrixXd qWq_product = q * interactions * q.transpose();
  Eigen::MatrixXd weighted_product = qWq_product.cwiseProduct(alpha_outer_product);
  Eigen::MatrixXd hessian = weighted_product + weighted_product.transpose();
  return hessian;
}
