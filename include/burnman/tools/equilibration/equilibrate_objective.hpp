/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_OBJ_HPP_INCLUDED
#define BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_OBJ_HPP_INCLUDED

#include <Eigen/Dense>
#include "burnman/tools/equilibration/equality_constraint_base.hpp"
#include "burnman/core/assemblage.hpp"

namespace burnman {
namespace equilibration {

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
 * @param x Parameter vector.
 * @param assemblage Taget assemblage (modified).
 * @param equality_constraints ConstraintGroup (vector of equality constraints).
 * @param reduced_composition_vector Amounts of independent elements.
 * @param reduced_free_composition_vector Amounts of independent elements in the free compositional vectors.
 * @return F(x) vector.
 */
Eigen::VectorXd F(
  const Eigen::VectorXd& x,
  Assemblage& assemblage,
  const ConstraintGroup& equality_constraints,
  const Eigen::VectorXd& reduced_composition_vector,
  const Eigen::MatrixXd& reduced_free_composition_vectors);

/**
 * @brief Jacobian function dF/dx for equilibration.
 *
 * @param x Parameter vector.
 * @param assemblage Target assemblage.
 * @param equality_constraints ConstraintGroup (vector of equality constraints).
 * @param reduced_free_compositional_vectors Reduced free compositional vectors.
 * @return Jacobian matrix, J(x).
 */
Eigen::MatrixXd J(
  const Eigen::VectorXd& x,
  Assemblage& assemblage, // TODO: could be const ref?
  const ConstraintGroup& equality_constraints,
  const Eigen::MatrixXd& reduced_free_composition_vectors);

} // namespace equilibration
} // namespace burnman

#endif // BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_OBJ_HPP_INCLUDED
