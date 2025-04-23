/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_SOLUTION_MODEL_HPP_INCLUDED
#define BURNMAN_CORE_SOLUTION_MODEL_HPP_INCLUDED

#include <string>
#include <map>
#include <Eigen/Dense>

/**
 * Base class for solution models.

  TODO
  
 */ 
class SolutionModel {

 public:

  virtual ~SolutionModel() = default;
  void process_solution_chemistry();

  // TODO:
  //  Re-order
  //  Check members and move to protected/private if possible...
  //  Consider size_t for n_endmembers etc.
  // Solution formulae -- vector with entry for each endmember
  std::vector<std::map<std::string, double>> solution_formulae;
  // Sites, list of lists with site possible atoms
  std::vector<std::vector<std::string>> sites;
  // Nested vector of occupancies for parsing
  std::vector<std::vector<std::vector<double>>> list_occupancies;
  int n_occupancies = 0;
  int n_endmembers;
  // endmembers --> vector of pointers to Mineral classes
  std::vector<std::string> formulas;
  // std::string empty_formula;
  // std::string general_formula;
  int n_sites;
  Eigen::ArrayXXd site_multiplicities;
  Eigen::ArrayXXd endmember_occupancies;
  Eigen::ArrayXXd endmember_n_occupancies;
  Eigen::ArrayXXd formula_multiplicities;

 protected:
  ;
 private:
  ;
};

#endif // BURNMAN_CORE_SOLUTION_MODEL_HPP_INCLUDED

// empty_formula [string] --> abbreviated formula with empty [] for sites
// general_formula --> formula with comma separated possible species
// n_sites [int]
// sites --> list of list of strings (check best way to implement, may ignore if unused)
// site_names --> vector of strings
// n_occupancies [int] --> sum of number of possible species
//  e.g. [[A][B], [B][C1/2D1/2]] would = 5 (2 on 1, 3 on 2)
// site_multiplicities [2D array] --> 1D for each endmember
// endmember_occupancies [2D array] --> as above
// endmember_noccupancies [2D array] --> number of atoms instead of fraction

// TODO:
//  Base class, ideal, asymmetric
//  Public getters and protected compute functions
//  Make virtual where needed to override in derived classes
/*
excess_gibbs_free_energy
excess_partial_gibbs_free_energies
excess_volume
excess_partial_volumes
excess_entropy
excess_partial_entropies
excess_enthalpy
Cp_excess
alphaV_excess
VoverKT_excess
// In derived
activity_coefficients
activities
// Asymmetric
Constructor logic
phi
non_ideal_interactions
non_ideal_excess_partial_gibbs
excess_partial_gibbs_free_energies --> override
excess_partial_entropies --> override
excess_partial_volumes --> override
gibbs_hessian
entorpy_hessian
volume_hessian
*/