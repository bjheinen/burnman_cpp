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
#include <iostream>

#include "burnman/core/solution_model.hpp"
#include "burnman/utils/constants.hpp"
#include "burnman/utils/string_utils.hpp"

void SolutionModel::process_solution_chemistry() {


  // Set class n_endmembers and n_sites
  n_endmembers = formulas.size();
  n_sites = std::count(formulas[0].begin(), formulas[0].end(), '[');

  // Check that number of sites is constant
  for (const std::string& f : formulas) {
    if (std::count(f.begin(), f.end(), '[') != n_sites) {
      throw std::runtime_error("All formulae must have the same number of distinct sites.");
    }
  }

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

  // TODO: move these to more senible place, maybe factor our into utils
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
  endmember_noccupancies = endmember_occupancies * site_multiplicities;

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
  std::string empty_formula;
  std::string general_formula;
  std::regex_replace(formulas[0], std::regex(R"(\[.*?\])"), "[]");
  // Split empty formula on [, then re-join with site species for general
  std::regex split_brackets(R"(\[)");
  std::sregex_token_iterator it_empty_split(
    empty_formula.begin(),
    empty_formula.end(),
    split_brackets, -1);
  std::vector<std::string> split_empty(++it_empty_split, {});
  general_formula = split_empty[0];
  for (int i = 0; i < n_sites; ++i) {
    general_formula += "[" + utils::join(sites[i], ",") + split_empty[i + 1];
  }

  // solution_model attributes to check:
  /*
    site_names
    solution_formulae
    n_sites
    sites
    site_multiplicities
    n_occupancies
    endmember_occupancies
    endmember_noccupancies
    empty_formula
    general_formula
  */

}
