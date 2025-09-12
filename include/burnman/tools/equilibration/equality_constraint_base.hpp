/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_TOOLS_EQUILIBRATION_EQUALITY_CONSTRAINT_BASE_HPP_INCLUDED
#define BURNMAN_TOOLS_EQUILIBRATION_EQUALITY_CONSTRAINT_BASE_HPP_INCLUDED

#include <memory>
#include <vector>
#include <Eigen/Dense>

namespace burnman {
namespace equilibration {

// Forward declaration
class Assemblage;

/**
 * @brief Base class for linear equality constraints.
 *
 * Use `make_constraint<ConstraintType>()' to construct a constraint.
 * Currently implemented constraints are:
 *   PressureConstraint
 *   TemperatureConstraint
 *   EntropyConstraint
 *   VolumeConstraint
 *   PTEllipseConstraint
 *   LinearXConstraint
 *   PhaseFractionConstraint
 *   PhaseCompositionConstraint
 *
 * Use `EqualityConstraint::evaluate()' to compute F.
 * Use `EqualityConstraint::derivative()' to compute J.
 *
 * @note Equality constraints should implement evaluate() and
 * derivative() functions, as well as a clone() method.
 */
class EqualityConstraint {
public:
  virtual ~EqualityConstraint() = default;
  virtual std::unique_ptr<EqualityConstraint> clone() const = 0;
  virtual double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const = 0;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const = 0;
};

/**
 * @brief Grouped top-level constraints.
 *
 * Use ConstraintGroup to keep expanded constraint vectors together.
 * Single constraints should also be put in a ConstraintGroup for
 * normalisation.
 *
 * @see `make_constraints_from_array' and `wrap_constraint'.
 */
using ConstraintGroup = std::vector<std::unique_ptr<EqualityConstraint>>;

/**
 * @brief Nested list of ConstraintGroups for equilibrate function.
 *
 * A ConstraintList should be constructed to pass to the equilibrate
 * function. Each top level ConstraintGroup can contain 1 or many
 * constraints. The equilibration routine will loop over all possible
 * lists of top level constraints constructed from the sub-constraints
 * in each group.
 *
 * @see `make_constraint_list'.
 */
using ConstraintList = std::vector<ConstraintGroup>;

} // namespace equilibration
} // namespace burnman

#endif // BURNMAN_TOOLS_EQUILIBRATION_EQUALITY_CONSTRAINT_BASE_HPP_INCLUDED
