/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/core/solution_models/solution_model_base.hpp"
#include <cstddef>
#include <algorithm>
#include <iterator>
#include <regex>
#include <stdexcept>
#include "burnman/utils/string_utils.hpp"

namespace burnman::solution_models {

SolutionModel::SolutionModel(const types::PairedEndmemberList& endmember_list) {
  this->n_endmembers = static_cast<Eigen::Index>(endmember_list.size());
  this->endmembers.reserve(static_cast<std::size_t>(this->n_endmembers));
  this->formulas.reserve(static_cast<std::size_t>(this->n_endmembers));
  // Unpack and store endmembers/formulas separately
  for (const auto& [mineral, formula] : endmember_list) {
    this->endmembers.push_back(mineral);
    this->formulas.push_back(formula);
  }
  // Process chemistry sets up site multiplicities/occupancies etc.
  process_solution_chemistry();
}

void SolutionModel::process_solution_chemistry() {

  // Reset n_occupancies count to zero
  this->n_occupancies = 0;

  // Set class n_sites
  this->n_sites = static_cast<Eigen::Index>(
    std::count(this->formulas[0].begin(), this->formulas[0].end(), '[')
  );
  // Check that number of sites is constant
  for (const std::string& f : this->formulas) {
    if (static_cast<Eigen::Index>(std::count(f.begin(), f.end(), '[')) != this->n_sites) {
      throw std::runtime_error("All formulae must have the same number of distinct sites.");
    }
  }

  // Make intermediate occ/mult data containers for parsing and resize as needed
  Eigen::ArrayXXd formula_multiplicities(this->n_endmembers, this->n_sites);
  std::vector<std::vector<std::vector<double>>> list_occupancies(
    static_cast<std::size_t>(this->n_endmembers)
  );
  this->solution_formulae.resize(static_cast<std::size_t>(this->n_endmembers));
  this->sites.resize(static_cast<std::size_t>(this->n_sites));
  // List occupancies is triple nested vector
  // std::vector<std::vector<std::vector<double>>>
  // Outer --> n_endmembers
  // Inner --> n_site
  // Innermost --> n_species (not known until parsed)
  for (std::size_t i = 0; i < static_cast<std::size_t>(this->n_endmembers); ++i) {
    list_occupancies[i].resize(static_cast<std::size_t>(this->n_sites));
    for (std::size_t j = 0; j < static_cast<std::size_t>(this->n_sites); ++j) {
      list_occupancies[i][j].resize(0);
    }
  }

  // TODO: move these to more sensible place, maybe factor out into string_utils
  // Regex strings for splits
  std::regex site_pattern_regex(R"(\[.*?\])");
  std::regex site_split_regex(R"(\[)");
  std::regex occ_split_regex(R"(\])");
  std::regex species_split_regex("[A-Z][^A-Z]*");
  std::regex species_frac_split_regex("([0-9][^A-Z]*)");

  // Loop over endmembers
  for (std::size_t i_mbr = 0; i_mbr < static_cast<std::size_t>(this->n_endmembers); ++i_mbr) {
    // Split formula into sites - 'Mg]3', 'Al]2', 'etc.'
    std::sregex_token_iterator it_sites(
      this->formulas[i_mbr].begin(), this->formulas[i_mbr].end(),
      site_split_regex, -1);
    // Discard string before first [ by ++it first
    std::vector<std::string> site_formulas(++it_sites, {});
    // Loop over sites in formula
    for (std::size_t i_site = 0; i_site < static_cast<std::size_t>(this->n_sites); ++i_site) {
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
      formula_multiplicities(
        static_cast<Eigen::Index>(i_mbr),
        static_cast<Eigen::Index>(i_site)
      ) = site_mult;

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
        this->solution_formulae[i_mbr][species_name] += site_mult * proportion_species_on_site;

        // Add to sites[i_site] if not present
        auto& site_species = this->sites[i_site];
        auto species_pos = std::find(site_species.begin(), site_species.end(), species_name);
        std::size_t i_el;
        if (species_pos == site_species.end()) {
          // Use current size as index of next entry
          i_el = site_species.size();
          site_species.push_back(species_name);
          ++this->n_occupancies;
          // Append 0 to list_occupancies for already parsed endmembers
          // when found new species
          for (std::size_t k = 0; k < static_cast<std::size_t>(this->n_endmembers); ++k) {
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
  this->endmember_occupancies.resize(this->n_endmembers, this->n_occupancies);
  this->site_multiplicities.resize(this->n_endmembers, this->n_occupancies);

  // Loop over endmembers again
  for (std::size_t i_mbr = 0; i_mbr < static_cast<std::size_t>(this->n_endmembers); ++i_mbr) {
    Eigen::Index n_species = 0;
    for (std::size_t i_site = 0; i_site < static_cast<std::size_t>(this->n_sites); ++i_site) {
      for (std::size_t i_el = 0; i_el < list_occupancies[i_mbr][i_site].size(); ++i_el) {
        this->endmember_occupancies(
          static_cast<Eigen::Index>(i_mbr),
          n_species
        ) = list_occupancies[i_mbr][i_site][i_el];
        this->site_multiplicities(
          static_cast<Eigen::Index>(i_mbr),
          n_species
        ) = formula_multiplicities(
          static_cast<Eigen::Index>(i_mbr),
          static_cast<Eigen::Index>(i_site)
        );
        ++n_species;
      }
    }
  } // Endmember loop

  // Element-wise multiply for n_occ
  this->endmember_n_occupancies = this->endmember_occupancies * this->site_multiplicities;

  // Get site names
  this->site_names.clear();
  for (std::size_t i_site = 0; i_site < static_cast<std::size_t>(this->n_sites); ++i_site) {
    // Grab uppercase letter from index (works to 26!)
    char site_id = static_cast<char>('A' + i_site);
    for (const std::string& sp : this->sites[i_site]) {
      this->site_names.push_back(sp + "_" + site_id);
    }
  }

  // NOTE:
  //   Currently A[Mg]2SiO4, [Fe]2SiO4 --> A[Mg,Fe]2SiO4, but
  //   [Mg]2SiO4, A[Fe]2SiO4 --> [Mg,Fe]2SiO4 -- does this matter?
  // Do parsed chemical formula etc. (general and empty)
  // Replace [sites] with [] on one formula for empty formula
  this->empty_formula = std::regex_replace(this->formulas[0], site_pattern_regex, "[]");
  // Split original formula on [.*?] site pattern
  std::sregex_token_iterator it_empty_split(
    this->formulas[0].begin(),
    this->formulas[0].end(),
    site_pattern_regex, -1);
  std::sregex_token_iterator end;
  std::vector<std::string> split_empty(it_empty_split, end);
  // Take first token to start general formula (may be empty)
  this->general_formula = split_empty[0];
  // Loop through and replace sites with multi species lists
  for (std::size_t i = 0; i < static_cast<std::size_t>(this->n_sites); ++i) {
    this->general_formula += "[" + utils::join(this->sites[i], ",") + "]";
    // Append final part after site if present
    if (i + 1 < split_empty.size()) {
      this->general_formula += split_empty[i + 1];
    }
  }
}

Eigen::Index SolutionModel::get_n_endmembers() const {
  return this->n_endmembers;
}

Eigen::Index SolutionModel::get_n_sites() const {
  return this->n_sites;
}

Eigen::Index SolutionModel::get_n_occupancies() const {
  return this->n_occupancies;
}

const Eigen::ArrayXXd& SolutionModel::get_site_multiplicities() const {
  return this->site_multiplicities;
}

const Eigen::ArrayXXd& SolutionModel::get_endmember_occupancies() const {
  return this->endmember_occupancies;
}

const Eigen::ArrayXXd& SolutionModel::get_endmember_n_occupancies() const {
  return this->endmember_n_occupancies;
}

const std::vector<std::string>& SolutionModel::get_site_names() const {
  return this->site_names;
}

const std::string& SolutionModel::get_empty_formula() const {
  return this->empty_formula;
}

const std::string& SolutionModel::get_general_formula() const {
  return this->general_formula;
}

const std::vector<std::string>& SolutionModel::get_formulas() const {
  return this->formulas;
}

const std::vector<std::vector<std::string>>& SolutionModel::get_sites() const {
  return this->sites;
}

const std::vector<std::map<std::string, double>>& SolutionModel::get_solution_formulae() const {
  return this->solution_formulae;
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

} // namespace burnman::solution_models
