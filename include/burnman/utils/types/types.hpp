/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_TYPES_HPP_INCLUDED
#define BURNMAN_UTILS_TYPES_HPP_INCLUDED

#include "burnman/utils/types/simple_types.hpp"
#include "burnman/utils/types/ndarray.hpp"
#include "burnman/utils/types/mineral_params.hpp"
#include "burnman/eos/components/excess_params.hpp"

namespace burnman {
  /**
  * @namespace burnman::types
  * @brief Custom types used in burnman
  *
  * TODO: docs on main user types - FormulaMap, MineralParams,
  * ExcessParams, etc. Also NDArray etc.
  */
  namespace types {
    // Aliases to expose excess_params types to user via burnman::types
    using LandauParams = burnman::eos::excesses::LandauParams;
    using LandauSLB2022Params = burnman::eos::excesses::LandauSLB2022Params;
    using LandauHPParams = burnman::eos::excesses::LandauHPParams;
    using LinearParams = burnman::eos::excesses::LinearParams;
    using BraggWilliamsParams = burnman::eos::excesses::BraggWilliamsParams;
    using MagneticChsParams = burnman::eos::excesses::MagneticChsParams;
    using DebyeParams = burnman::eos::excesses::DebyeParams;
    using DebyeDeltaParams = burnman::eos::excesses::DebyeDeltaParams;
    using EinsteinParams = burnman::eos::excesses::EinsteinParams;
    using EinsteinDeltaParams = burnman::eos::excesses::EinsteinDeltaParams;
    using ExcessParamVector = burnman::eos::excesses::ExcessParamVector;
  } // namespace types
} // namespace burnman

#endif // BURNMAN_UTILS_TYPES_HPP_INCLUDED
