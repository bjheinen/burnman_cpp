/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_IMPL_HPP_INCLUDED
#define BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_IMPL_HPP_INCLUDED

#include <vector>
#include <Eigen/Dense>
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/tools/equilibration/equality_constraint_base.hpp"
#include "burnman/tools/equilibration/equilibrate_types.hpp"
#include "burnman/core/assemblage.hpp"

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
 * @param composition Bulk composition.
 * @param assemblage Target assemblage to equilibrate.
 * @param equality_constraints Vector of top level ConstraintGroups.
 * @param free_compositional_vectors Optional element maps for bulk compositional degrees of freedom.
 * @param tol Convergence tolerance.
 * @param store_iterates Store parameter vector and function value convergence history.
 * @param store_assemblage Store assemblage object copies in result.
 * @param max_iterations Maximum allowed iterations.
 * @param verbose Print extra iteration info.
 * @return EquilibrateResult with solution and parameters.
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

#endif // BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_IMPL_HPP_INCLUDED
