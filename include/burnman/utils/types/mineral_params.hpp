/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_TYPES_MINERAL_PARAMS_HPP_INCLUDED
#define BURNMAN_UTILS_TYPES_MINERAL_PARAMS_HPP_INCLUDED

#include <optional>
#include <string>
#include "burnman/utils/types/simple_types.hpp"

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
  std::optional<double> T_einstein; // Used for SLB --> can be calculated from S_0 and napfu

  // Electronic parameters (SLB Conductive)
  std::optional<double> bel_0;
  std::optional<double> gel;

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

#endif // BURNMAN_UTILS_TYPES_MINERAL_PARAMS_HPP_INCLUDED
