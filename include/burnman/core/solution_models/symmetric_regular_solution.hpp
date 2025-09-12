/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_SOLUTION_MODELS_SYMMETRIC_REGULAR_HPP_INCLUDED
#define BURNMAN_CORE_SOLUTION_MODELS_SYMMETRIC_REGULAR_HPP_INCLUDED

#include <vector>
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/core/solution_models/asymmetric_regular_solution.hpp"

namespace burnman {
namespace solution_models {

/**
 * @class SymmetricRegularSolution
 * @brief Convenience class for a symmetric regular solution model.
 *
 * This class is a special case of `AsymmetricRegularSolution' with all
 * alphas set to 1.
 *
 */
class SymmetricRegularSolution : public AsymmetricRegularSolution{
 public:
  SymmetricRegularSolution(
    const types::PairedEndmemberList& endmember_list,
    std::vector<std::vector<double>> energy_interaction,
    std::vector<std::vector<double>> volume_interaction = {},
    std::vector<std::vector<double>> entropy_interaction = {});
};

} // namespace solution_models
} // namespace burnman

#endif // BURNMAN_CORE_SOLUTION_MODELS_SYMMETRIC_REGULAR_HPP_INCLUDED
