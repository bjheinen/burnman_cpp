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

#include <string>
#include <utility>
#include <unordered_map>
#include <vector>
#include <Eigen/Dense>
#include "burnman/core/assemblage.hpp"
#include "burnman/utils/types.hpp"

// namespace equilibrate maybe?

// TODO: move to types.hpp?
using FreeVectorMap = std::unordered_map<std::string, double>;

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
 * @brief Holds result of equilibration along with parameters
 */
struct EquilibrateResult {
  NDArray<DampedNewtonResult> sol_array;
  EquilibrationParameters prm;
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

/**
 * @brief Builds the equilibration parameter object.
 */
EquilibrationParameters get_equilibration_parameters(
  const Assemblage& assemblage,
  const FormulaMap& composition,
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
 */
std::pair<Eigen::MatrixXd, Eigen::VectorXd> calculate_constraints(
  const Assemblage& assemblage,
  int n_free_compositional_vectors);

/**
 * @brief Objective function for equilibration.
 *
 * The vector-valued function for which the root is sought.
 * The first two vector values depend on the equality_constraints chosen:
 *   eq[i] is PressureConstraint, F[i] = P - eq[i].value
 *   eq[i] is TemperatureConstraint, F[i] = T - eq[i].value
 *   eq[i] is EntropyConstraint, F[i] = entropy - eq.value
 *   eq[i] is VolumeConstraint, F[i] = volume - eq.value
 *   eq[i] is PTEllipseConstraint, F[i] = norm(([P, T] - eq.centre)/eq.scaling) - 1
 *   eq[i] is LinearXConstraint, F[i] = eq.A.dot(x) - eq.b
 * See `EqualityConstraint::evaluate()' for details.
 *
 * The next set of vector values correspond to the reaction affinities.
 * The final set of vector values correspond to the bulk composition constraints.
 *
 * TODO: @param
 */
Eigen::VectorXd F(
  const Eigen::VectorXd& x,
  Assemblage& assemblage,
  const std::vector<std::unique_ptr<EqualityConstraint>>& equality_constraints,
  const Eigen::VectorXd& reduced_composition_vector,
  const Eigen::MatrixXd& reduced_free_composition_vectors);

/**
 * @brief Jacobian function for equilibration, dF/dx.
 *
 * TODO: @param
 *
 */
Eigen::MatrixXd J(
  const Eigen::VectorXd& x,
  Assemblage& assemblage,
  const std::vector<std::unique_ptr<EqualityConstraint>>& equality_constraints,
  const Eigen::MatrixXd& reduced_free_composition_vectors);

/**
 * @brief Finds equilibrium state of assemblage subject to constraints.
 *
 * This finds the thermodynamic equilibrium state of an assemblage subject to
 * given equality constraints by solving a set of nonlinear equations
 * related to the chemical potentials and other state variables of the system.
 *
 * Usage requires an assemblage and 2 + n_c equality constraints, where n_c is
 * the number of bulk compositional degrees of freedom. The equilibrate
 * function attempts to find the remaining unknowns that satisfy those constraints.
 *
 * See `EqualityConstraint' for information on constraints TODO: note namespace here
 * Possible constraints are:
 *   PressureConstraint
 *   TemperatureConstraint
 *   EntropyConstraint
 *   VolumeConstraint
 *   PTEllipseConstraint
 *   PhaseFractionConstraint
 *   PhaseCompositionConstraint
 *   LinearXConstraint
 *
 * TODO: @param
 */
EquilibrateResult equilibrate(
  const FormulaMap& composition,
  Assemblage& assemblage,
  const ConstraintList& equality_constraints,
  const std::vector<FreeVectorMap>& free_compositional_vectors = {},
  double tol = 1.0e-3,
  bool store_iterates = false,
  bool store_assemblage = true,
  int max_iterations = 100,
  bool verbose = false);

#endif // BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_HPP_INCLUDED
