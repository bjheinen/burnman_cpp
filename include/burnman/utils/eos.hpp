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

#include <string>
#include <unordered_map>
#include <optional>
#include <variant>
#include <vector>

// Type alias for formula map
using FormulaMap = std::unordered_map<std::string, int>;

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
    int n, factor;
    double Wh, Wv, deltaH, deltaV;
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
 * Enum used to define EOS Type
 */
enum class EOSType {
  Auto, // Used to set EOS from params
  Custom, // Used when user passes custom EOSType
  Vinet,
  BM3,
  BM2,
  MGD2,
  MGD3
};

/**
 * Four Cp parameters used in HP and CORK
 */
struct CpParams {
  double a, b, c, d;
};

/**
 * Parameters for CORK
 */
struct CorkParams {
  double a_0, a_1, b, c_0, c_1, d_0, d_1;
};

/**
 * Struct to hold mineral EOS parameters
 */
struct MineralParams {
  // Using std::optional so one struct can be used for all EOS types
  // We want to be able to switch between EOS without redefining minerals
  // Sometime we want to calculate default values (e.g. P0, T0) if not set
  // Should remove std::optional for anything always required (V0?, etc.)

  // Mineral properties
  std::optional<std::string> name;
  std::optional<FormulaMap> formula;
  std::optional<EOSType> equation_of_state;
  std::optional<int> napfu;
  std::optional<double> molar_mass;

  // Reference conditions
  std::optional<double> P_0;
  std::optional<double> T_0;
  std::optional<double> E_0;
  std::optional<double> F_0;

  // Reference volume
  std::optional<double> V_0;

  // Elastic moduli
  std::optional<double> K_0;
  std::optional<double> Kprime_0;
  std::optional<double> Kdprime_0;
  std::optional<double> G_0;
  std::optional<double> Gprime_0;

  // Thermodynamic limit (infinite pressure) parameters
  std::optional<double> Kprime_inf;
  std::optional<double> Gprime_inf;

  // Other thermodynamic properties
  std::optional<double> S_0;
  std::optional<double> H_0;
  std::optional<double> Cv; // DKS
  std::optional<CpParams> Cp;

  // Thermal parameters
  std::optional<double> debye_0;
  std::optional<double> grueneisen_0;
  std::optional<double> q_0;
  std::optional<double> a_0; // alpha_0
  std::optional<double> dKdT_0;
  std::optional<double> m; // DKS free param
  std::optional<double> a; // DKS free param
  std::optional<double> eta_s_0; // for shear strain derivative of grueneisen parameter (SLB, DKS)

  // Order of expansion (DKS)
  std::optional<int> order_theta;
  std::optional<int> order_f;

  // DKS free parameters for vol dependence of Tel and zeta
  // Tel - lower conductivity limit
  // zeta -  electronic heat capacity coefficient
  std::optional<double> Tel_0;
  std::optional<double> zeta_0;
  std::optional<double> xi;
  std::optional<double> eta;

  // Brosh-CALPHAD fitting parameters
  std::optional<double> delta_0;
  std::optional<double> delta_1;
  std::optional<double> b_0;
  std::optional<double> b_1;

  // CORK params
  std::optional<double> cork_T;
  std::optional<double> cork_P;
  std::optional<CorkParams> cork_params;
};

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
   * Struct for GSL Brent root finding
   * Used in Bragg Williams excess function to find
   * Q in Gibbs calculation.
   */
  struct BWReactParams {
    double delta_H, temperature, W;
    int n, f_0, f_1;
  };

} // End namespace ParamGSL

#endif // BURNMAN_UTILS_EOS_HPP_INCLUDED
