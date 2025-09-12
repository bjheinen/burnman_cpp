/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_UTILS_HPP_INCLUDED
#define BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_UTILS_HPP_INCLUDED

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/tools/equilibration/equilibrate_types.hpp"
#include "burnman/core/assemblage.hpp"

namespace burnman {
namespace equilibration {

/**
 * @brief Builds the equilibration parameter object.
 *
 * @param assemblage Target assemblage.
 * @param composition Bulk composition of the assemblage
 * @param free_compositional_vectors Optional free compositional vectors.
 * @return EquilibrationParameters object.
 */
EquilibrationParameters get_equilibration_parameters(
  const Assemblage& assemblage,
  const types::FormulaMap& composition,
  const std::vector<std::unordered_map<std::string, double>>& free_compositional_vectors);

/**
 * @brief Calculates the linear inequality constraints bounding the valid parameter space for an assemblage.
 *
 * The constraints are:
 *   - Pressure and temperature must be +ve.
 *   - All phase fractions must be +ve.
 *   - All site-species occupancies must be +ve.
 *
 * The constraints are stored in a vector (b) and matrix (A).
 * The sign convention is chosen such that the constraint is satisfied
 * if A·x + b < eps.
 *
 * @param assemblage Target assemblage.
 * @param n_free_compositional_vectors Number of bulk compositional degrees of freedom.
 */
std::pair<Eigen::MatrixXd, Eigen::VectorXd> calculate_constraints(
  const Assemblage& assemblage,
  int n_free_compositional_vectors);

/**
 * @brief Makes the starting parameter vector for equilibrium problem.
 *
 * The parameters are:
 *   - pressure
 *   - temperature
 *   - absolute amount of each phase
 *       if a phase if a solution with > 1 endmember, the following parameters
 *       are the mole fractions of the independent endmembers in the solution
 *       (i.e. missing the first endmember as they must sum to one).
 *
 * @param assemblage The target Assemblage object
 * @param n_free_compositional_vectors Number of bulk compositional degrees of freedom.
 * @return Parameter vector
 */
Eigen::VectorXd get_parameter_vector(
  const Assemblage& assemblage,
  int n_free_compositional_vectors = 0);

/**
 * @brief Returns the absolute amounts of all the endmembers in the solution.
 */
Eigen::VectorXd get_endmember_amounts(
  const Assemblage& assemblage);

/**
 * @brief Sets state and composition of an Assemblage from parameter vector.
 *
 * Parameter vector contains P, T, X.
 */
void set_composition_and_state_from_parameters(
  Assemblage& assemblage,
  const Eigen::VectorXd& parameters);

} // namespace equilibration
} // namespace burnman

#endif // BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_UTILS_HPP_INCLUDED
