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
 * @brief Contains parameters and data required for equilibration.
 */
struct EquilibrationParameters {
  std::vector<std::string> parameter_names;           ///< Names of all parameters.
  Eigen::VectorXd bulk_composition_vector;            ///< Bulk composition vector of assemblage elements.
  Eigen::VectorXd reduced_composition_vector;         ///< Bulk composition resitricted to independent elements.
  Eigen::MatrixXd free_compositional_vectors;         ///< Matrix of free compositional vectors (each row one vector).
  Eigen::MatrixXd reduced_free_compositional_vectors; ///< Free compositional vectors restricted to independent elements.
  Eigen::VectorXd constraint_vector;                  ///< Constraint vector, b, for constraints: (A·x + b).
  Eigen::MatrixXd constraint_matrix;                  ///< Constraint matrix, A, for contraints (A·x + b).
  Eigen::ArrayXi phase_amount_indices;                ///< Indices in parameter vector of phase amounts.
  Eigen::Index n_parameters;                          ///< Number of parameters for the equilibrium problem.
};

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
