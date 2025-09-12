/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_TYPES_GSL_PARAMS_HPP_INCLUDED
#define BURNMAN_UTILS_TYPES_GSL_PARAMS_HPP_INCLUDED

#include "burnman/utils/types/mineral_params.hpp"

// TODO: move to /eos/components - no need to expose as not user-facing
namespace burnman {
namespace eos {

/**
 * @namespace burnman::eos::gsl_params
 * @brief Structs to hols parameters when making GSL function objects
 */
namespace gsl_params {
  /**
   * Struct for GSL Brent root finding
   * Used for volume finding in EOS where only pressure needed
   * as an additional argument:
   *   bm, bm4, vinet, macaw, morse_potential, spock
   */
  struct SolverParams_P {
    const MineralParams& params;
    double pressure;
  };
  /**
   * Struct for GSL Brent root finding
   * Used for volume finding in EOS where pressure and
   * temperature needed as additional arguments:
   *   mgd3
   */
  struct SolverParams_PT {
    const MineralParams& params;
    double pressure;
    double temperature;
  };
  /**
   * Struct for GSL root finding
   * Used for volume finding in EOS where
   * P, T needed along with SLB specific params.
   *  SLB2, SLB3, etc.
   */
  struct SolverParams_SLB {
    const MineralParams& params;
    double pressure;
    double temperature;
    double a1_ii, a2_iikk;
    double b_iikk, b_iikkmm;
    double bel_0, gel;
  };
  /**
   * Struct for GSL Brent root finding
   * Used in Bragg Williams excess function to find
   * Q in Gibbs calculation.
   */
  struct BWReactParams {
    double delta_H, temperature, W;
    int n;
    double f_0, f_1;
  };

} // namespace gsl_params
} // namespace eos
} // namespace burnman

#endif // BURNMAN_UTILS_TYPES_GSL_PARAMS_HPP_INCLUDED
