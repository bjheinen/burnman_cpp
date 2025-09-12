/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_EOS_TYPES_EXCESS_PARAMS_HPP_INCLUDED
#define BURNMAN_EOS_TYPES_EXCESS_PARAMS_HPP_INCLUDED

#include <variant>
#include <vector>

namespace burnman {
namespace eos {
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

} // namespace excesses
} // namespace eos
} // namesapce burnman

#endif // BURNMAN_EOS_TYPES_EXCESS_PARAMS_HPP_INCLUDED
