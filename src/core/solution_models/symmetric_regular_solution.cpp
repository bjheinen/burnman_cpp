/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/core/solution_models/symmetric_regular_solution.hpp"

namespace burnman::solution_models {

// SymmetricRegularSolution Constructor
SymmetricRegularSolution::SymmetricRegularSolution(
  const PairedEndmemberList& endmember_list,
  std::vector<std::vector<double>> energy_interaction,
  std::vector<std::vector<double>> volume_interaction,
  std::vector<std::vector<double>> entropy_interaction)
: AsymmetricRegularSolution(
    endmember_list,
    std::vector<double>(endmember_list.size(), 1.0),
    energy_interaction,
    volume_interaction,
    entropy_interaction) {}

} // namespace burnman::solution_models
