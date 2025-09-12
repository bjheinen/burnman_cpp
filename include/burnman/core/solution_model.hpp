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

// Utility header to include all solution models
#include "burnman/core/solution_models/solution_model_base.hpp"
#include "burnman/core/solution_models/ideal_solution.hpp"
#include "burnman/core/solution_models/asymmetric_regular_solution.hpp"
#include "burnman/core/solution_models/symmetric_regular_solution.hpp"

namespace burnman {
  /**
  * @namespace burnman::solution_models
  * @brief Contains classes and functions for solid solutions models.
  *
  * This namespace includes base and derived solution model classes.
  * All `burnman::Solution' mineral objects use a solution model to
  * define how the endmembers interact.
  * These classes handle chemical formulae, sites and site occupancies,
  * as well as excess thermodynamic properties.
  *
  * Currently implemented solution models are:
  * - IdealSolution
  * - AsymmetricRegularSolution
  * - SymmetricRegularSolution
  */
  namespace solution_models {
    // blank for documentation only
  } // namespace solution_models
} // namespace burnman

#endif // BURNMAN_CORE_SOLUTION_MODEL_HPP_INCLUDED
