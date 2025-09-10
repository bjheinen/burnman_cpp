/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_TOOLS_EQUILIBRATION_EQUALITY_CONSTRAINT_HELPERS_HPP_INCLUDED
#define BURNMAN_TOOLS_EQUILIBRATION_EQUALITY_CONSTRAINT_HELPERS_HPP_INCLUDED

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include "burnman/tools/equilibration/equality_constraint_base.hpp"
#include "burnman/tools/equilibration/equality_constraint_variants.hpp"

// TODO: Currently cant make a constraint group for LinearXConstraint,
//       PhaseFractionConstraint, or PhaseCompositionConstraint.
//       Can add extra template specialisation for these, but probably not needed ever.

/**
 * @brief Constructs a constraint object.
 *
 * This factory function generates a constraint object from
 * a set of data. See indivudal EqualityConstraint defintions for
 * required arguments.
 *
 * Usage:
 *   std::unique_ptr<EqualityConstraint> p_constraint =
 *     make_constraint<PressureConstraint>(25.0e9);
 *
 * @tparam ConstraintT Type of constraint to create (derived from EqualityConstraint).
 * @param args Arguments forwarded to the constructor of ConstraintT.
 * @return Unique pointer to the constraint object.
 */
template <typename ConstraintT, typename... Args>
std::unique_ptr<EqualityConstraint> make_constraint(Args&&... args) {
  return std::make_unique<ConstraintT>(std::forward<Args>(args)...);
}

/**
 * @brief Constructs a group of constraints from vector of scalar values.
 *
 * This factory function constructs a ConstraintGroup (a std::vector of
 * std::unique_ptr<EqualityConstraint>) by making one constraint per value
 * in the input vector.
 *
 * @tparam ConstraintT Type of constraint to create (derived from EqualityConstraint).
 * @param values Vector of scalar values used to construct individual constraints.
 * @return ConstraintGroup containing one constraint per value.
 */
template <typename ConstraintT>
ConstraintGroup make_constraints_from_array(
  const Eigen::ArrayXd& values
) {
  std::vector<std::unique_ptr<EqualityConstraint>> constraints;
  constraints.reserve(static_cast<std::size_t>(values.size()));
  for (double v : values) {
    constraints.push_back(make_constraint<ConstraintT>(v));
  }
  return constraints;
}

/**
 * @brief Specialisation for creating PTEllipseConstraint objects from centre/scale arrays.
 *
 * Constructs a ConstraintGroup of PTEllipseConstraint objects from a pair of 
 * Eigen::ArrayXd representing centres and scales of ellipses.
 *
 * @param values Pair of arrays: first = centres, second = scales
 * @return ConstraintGroup of PTEllipseConstraint objects
 */
template <>
ConstraintGroup make_constraints_from_array<PTEllipseConstraint>(
  const std::pair<Eigen::ArrayXd, Eigen::ArrayXd>& values
) {
  const Eigen::ArrayXd& centres = values.first;
  const Eigen::ArrayXd& scales = values.second;
  std::vector<std::unique_ptr<EqualityConstraint>> constraints;
  constraints.reserve(static_cast<std::size_t>(centres.size()));
  for (Eigen::Index i = 0; i < centres.szie(); ++i) {
    constraints.push_back(make_constraint<PTEllipseConstraint>(centres(i), scales(i)));
  }
  return constraints;
}


// Helpers to wrap single constraint or pass through an existing group
/**
 * @brief Wrap a single constraint into a ConstraintGroup.
 *
 * This is useful when constructing a ConstraintList from a mix of
 * single constraints and pre-existing groups. The single constraint
 * will be placed as the only element in a ConstraintGroup.
 *
 * @param c The unique pointer to a single EqualityConstraint.
 * @return ConstraintGroup containing the single constraint.
 */
inline ConstraintGroup wrap_constraint(std::unique_ptr<EqualityConstraint> c) {
  return {std::move(c)};
}

/**
 * @brief Pass through an existing ConstraintGroup.
 *
 * This overload allows wrap_constraint to accept a ConstraintGroup
 * without modification, which simplifies helpers that
 * construct a ConstraintList from mixed single/group inputs.
 *
 * @param g A pre-existing ConstraintGroup.
 * @return The same ConstraintGroup (unchanged).
 */
inline ConstraintGroup wrap_constraint(ConstraintGroup g) {
  return g;
}

/**
 * @brief Construct a ConstraintList from a sequence of constraints or constraint groups.
 *
 * This helper builds a nested ConstraintList from any mix of:
 *   - Single `EqualityConstraint` objects (as `std::unique_ptr<EqualityConstraint>`)
 *   - Pre-existing `ConstraintGroup`s
 * Each argument is automatically wrapped into a `ConstraintGroup` if needed.
 *
 * @param args The constraints or constraint groups to include in the list.
 * @return `ConstraintList` containing all the constraints/groups.
 *
 */
template <typename... Args>
ConstraintList make_constraint_list(Args&&... args) {
  ConstraintList list;
  (list.push_back(wrap_constraint(std::forward<Args>(args))), ...);
  return list;
}

#endif // BURNMAN_TOOLS_EQUILIBRATION_EQUALITY_CONSTRAINT_HELPERS_HPP_INCLUDED
