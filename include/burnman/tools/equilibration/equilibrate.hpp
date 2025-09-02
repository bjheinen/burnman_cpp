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

#include <utility>
#include <vector>
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

/**
 * @brief Custom function to return lambda bounds for damped Newton solver.
 *
 * Computes bounds on lambda for the damped affine invariant modification to
 * Newton's method for nonlinear problems (Deuflhard, 1974;1975;2004).
 *
 * @param dx The proposed Newton step.
 * @param x The parameter vector.
 *
 * @returns (min, max)
 */
std::pair<double, double> lambda_bounds_func(
  const Eigen::VectorXd& dx,
  const Eigen::VectorXd& x,
  const std::vector<int>& endmembers_per_phase);

#endif // BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_HPP_INCLUDED
