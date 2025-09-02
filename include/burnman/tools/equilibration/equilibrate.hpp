/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_HPP_INCLUDED
#define BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_HPP_INCLUDED

#include <Eigen/Dense>
#include "burnman/core/assemblage.hpp"

// namespace equilibrate maybe?

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
 * @param n_free_compositional_vectors
 * @return Parameter vector
 */
Eigen::VectorXd get_parameter_vector(
  const Assemblage& assemblage,
  int n_free_compositional_vectors = 0);

/**
 * @brief Get the absolute amounts of all the endmembers in the solution.
 */
Eigen::VectorXd get_endmember_amounts(
  const Assemblage& assemblage);

/**
 * @brief Helper function to set state and composition from parameter vector.
 *
 * Parameter vector contains P, T, X.
 */
void set_composition_and_state_from_parameters(
  Assemblage& assemblage,
  const Eigen::VectorXd& parameters);

#endif // BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_HPP_INCLUDED
