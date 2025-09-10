/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_EOS_HPP_INCLUDED
#define BURNMAN_UTILS_EOS_HPP_INCLUDED

#include <variant>
#include <vector>
#include "burnman/utils/types/simple_types.hpp"

/**
 * @namespace excesses
 * Used for calculating thermodynamic property modifiers
 */
namespace excesses {
  /**
  * Sruct to hold excesses data for EOS modifiers
  */
  struct Excesses {
      double G = 0.0;
      double dGdT = 0.0;
      double dGdP = 0.0;
      double d2GdT2 = 0.0;
      double d2GdP2 = 0.0;
      double d2GdPdT = 0.0;
      // Overloaded += operator to add Excesses
      Excesses& operator+=(const Excesses& other) {
          G += other.G;
          dGdT += other.dGdT;
          dGdP += other.dGdP;
          d2GdT2 += other.d2GdT2;
          d2GdP2 += other.d2GdP2;
          d2GdPdT += other.d2GdPdT;
          return *this;
      }
  };

  /**
  * Parameters for tricritical landau correction
  * following Putnis (1992)
  */
  struct LandauParams {
    double Tc_0, V_D, S_D;
  };

  /**
  * Parameters for tricritical landau correction
  * following Stixrude & Lithgow-Bertelloni (2022)
  */
  struct LandauSLB2022Params {
    double Tc_0, V_D, S_D;
  };

  /**
  * Parameters for tricritical landau correction
  * following Holland & Powell (1998)
  */
  struct LandauHPParams {
    double T_0, P_0, Tc_0, V_D, S_D;
  };

  /**
  * Parameters for Darken's quadratic correction
  * following Powell (1987)
  */
  struct LinearParams {
    double delta_V, delta_S, delta_E;
  };

  /**
  * Parameters for Bragg-Williams type correction (order-disorder)
  * following Holland & Powell (1987)
  */
  struct BraggWilliamsParams {
    int n;
    double factor, Wh, Wv, deltaH, deltaV;
  };

  /**
  * Parameters for magnetic contribution correction
  * following Chin, Hetzman & Sundman (1987) and Sundman (1991)
  */
  struct MagneticChsParams {
    double structural_parameter;
    double curie_T_0, curie_T_p;
    double magnetic_moment_0, magnetic_moment_p;
  };

  /**
  * Parameters excess contribution from Debye model
  */
  struct DebyeParams {
    double Cv_inf, Theta_0;
  };

  /**
  * Parameters for excess contribution from
  * thermal derivatives of Debye model
  */
  struct DebyeDeltaParams {
    double S_inf, Theta_0;
  };

  /**
  * Parameters excess contribution from Einstein model
  */
  struct EinsteinParams {
    double Cv_inf, Theta_0;
  };

  /**
  * Parameters for excess contribution from
  * thermal derivatives of Einstein model
  */
  struct EinsteinDeltaParams {
    double S_inf, Theta_0;
  };

  using ExcessParamVariant = std::variant<
    LandauParams,
    LandauSLB2022Params,
    LandauHPParams,
    LinearParams,
    BraggWilliamsParams,
    MagneticChsParams,
    DebyeParams,
    DebyeDeltaParams,
    EinsteinParams,
    EinsteinDeltaParams
  >;

  using ExcessParamVector = std::vector<ExcessParamVariant>;

} // End namespace ExcessParams

/**
 * Structs to hols parameters when making GSL function objects
 */
namespace ParamsGSL {
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

} // End namespace ParamGSL

#endif // BURNMAN_UTILS_EOS_HPP_INCLUDED
