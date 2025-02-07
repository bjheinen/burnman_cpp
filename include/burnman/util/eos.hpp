/*
  TODO: Copyright Notice!
*/
#ifndef BURNMAN_UTIL_EOS_HPP_INCLUDED
#define BURNMAN_UTIL_EOS_HPP_INCLUDED

#include <string>
#include <optional>

// Type alias for formula map
using FormulaMap = std::unordered_map<std::string, int>;

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
 * Four Cp parameters used in HP and CORK
 */
struct CpParams {
  double a, b, c, d;
}

/**
 * Parameters for CORK
 */
struct CorkParams {
  double a_0, a_1, b, c_0, c_1, d_0, d_1;
}


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
  std::optional<std::string> equation_of_state;
  std::optional<int> napfu;
  std::optional<double> molar_mass;

  // Reference conditions
  std::optional<double> P_0;
  std::optional<double> T_0;

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
  std::optional<double> eta_s_0; // for shear strain derivate of grueneisen parameter (SLB, DKS)

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
}


#endif // BURNMAN_UTIL_EOS_HPP_INCLUDED