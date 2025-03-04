/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <cmath>
#include "burnman/eos/property_modifiers.hpp"
#include "burnman/eos/debye.hpp"
#include "burnman/eos/einstein.hpp"
#include "burnman/optim/brent_solver.hpp"
#include "burnman/utils/constants.hpp"

// Anonymous namespace for BW helper functions
namespace {
  // Eq. A2-2
  double flnarxn(int n, double Q, int f_0, int f_1) {
    double log1mQ = std::log1p(-Q);
    return n / (n + 1)
      * (f_0 * (std::log(n) + log1mQ)
        + f_1 * log1mQ
        - f_0 * std::log1p(n * Q)
        - f_1 * std::log(n + Q));
  }

  double reaction_bragg_williams_gsl_wrapper(double Q, void* p) {
    // Cast pointer back to struct
    auto* params = static_cast<const ParamsGSL::BWReactParams*>(p);
    return params->delta_H
      + constants::physics::gas_constant * params->temperature
      * flnarxn(params->n, Q, params->f_0, params->f_1)
      + (2.0 * Q - 1.0) * params->W;
  }

  double order_gibbs(
    double pressure, double temperature,
    excesses::BraggWilliamsParams params,
    int f_0, int f_1
  ) {
    double W = params.Wh + pressure * params.Wv;
    double H_disord = params.deltaH + pressure * params.deltaV;

    // Make params struct for brent solver
    ParamsGSL::BWReactParams react_params{
      H_disord, temperature, W,
      params.n, f_0, f_1
    };
    // TODO: Check for non-converge in brent --> set Q=0 if error
    double Q = brent::find_root(
      &reaction_bragg_williams_gsl_wrapper,
      react_params,
      constants::precision::abs_tolerance,
      1.0 - constants::precision::abs_tolerance
    );
    // Reuse some identities
    double lognp1 = std::log1p(params.n);
    double logQm1 = std::log1p(-Q);
    double S = -constants::physics::gas_constant
      * (f_0
        * ((1 + params.n * Q) * (std::log1p(params.n*Q) - lognp1)
          + params.n * (1.0 - Q) * (std::log(params.n) + logQm1 - lognp1))
        + f_1
        * (params.n * (1.0 - Q) * (logQm1 - lognp1)
          + params.n * (params.n + Q) * (std::log(params.n + Q) - lognp1))
      ) / (params.n + 1);

    double G = (1.0 - Q) * H_disord + (1.0 - Q) * Q * W - temperature * S;
    // Ignoring Q for now...
    return G;
  }
}

namespace excesses {

  Excesses compute_excesses(
    double pressure,
    double temperature,
    LandauParams params
  ) {
    double Q2;
    double G;
    double dGdT;
    double dGdP;
    double d2GdT2;
    double d2GdP2;
    double d2GdPdT;
    double Tc = params.Tc_0 + params.V_D * pressure / params.S_D;
    double G_disordered = -params.S_D * ((temperature - Tc) + params.Tc_0 / 3.0);
    double dGdT_disordered = -params.S_D;
    double dGdP_disordered = params.V_D;
    if (temperature < Tc) {
      double Tc_0_over_Tc = params.Tc_0 / Tc;
      Q2 = std::sqrt(1.0 - temperature / Tc);
      G = params.S_D
        * ((temperature - Tc) * Q2 + params.Tc_0 * Q2 * Q2 * Q2 / 3.0)
        + G_disordered;
      dGdP = -params.V_D * Q2
        * (1.0 + 0.5 * temperature / Tc * (1.0 - Tc_0_over_Tc))
        + dGdP_disordered;
      dGdT = params.S_D * Q2
        * (1.5 - 0.5 * Tc_0_over_Tc)
        + dGdT_disordered;
      d2GdP2 = params.V_D * params.V_D * temperature
        / (params.S_D * Tc * Tc * Q2)
        * (temperature * (1.0 + Tc_0_over_Tc) / (4.0 * Tc)
          + Q2 * Q2 * (1.0 - Tc_0_over_Tc) - 1.0);
      d2GdT2 = -params.S_D / (Tc * Q2) * (0.75 - 0.25 * Tc_0_over_Tc);
      d2GdPdT = params.V_D / (2.0 * Tc * Q2)
        * (1.0 + (temperature / (2.0 * Tc) - Q2 * Q2) * (1.0 - Tc_0_over_Tc));
    } else {
      Q2 = 0.0;
      G = G_disordered;
      dGdT = dGdP_disordered;
      dGdP = dGdP_disordered;
      d2GdT2 = 0.0;
      d2GdP2 = 0.0;
      d2GdPdT = 0.0;
    }
    Excesses landau_ex{
      G,
      dGdT,
      dGdP,
      d2GdT2,
      d2GdP2,
      d2GdPdT};
    // double Q = std::sqrt(Q2);
    return landau_ex;
  }

  Excesses compute_excesses(
    double pressure,
    double temperature,
    LandauSLB2022Params params
  ) {
    double Q2;
    double G;
    double dGdT;
    double dGdP;
    double d2GdT2;
    double d2GdP2;
    double d2GdPdT;
    double Tc = params.Tc_0 + params.V_D * pressure / params.S_D;
    double G_disordered = -params.S_D * ((temperature - Tc) + params.Tc_0 / 3.0);
    double dGdT_disordered = -params.S_D;
    double dGdP_disordered = params.V_D;
    if (temperature < Tc) {
      Q2 = std::sqrt((Tc - temperature) / params.Tc_0);
      if (Q2 < 4.0) {
        G = params.S_D
          * ((temperature - Tc) * Q2 + params.Tc_0 * Q2 * Q2 * Q2 / 3.0)
          + G_disordered;
        dGdT = params.S_D * Q2 + dGdT_disordered;
        dGdP = -params.V_D * Q2 + dGdP_disordered;
        d2GdT2 = -0.5 * params.S_D / (Tc * Q2);
        d2GdP2 = -0.5 * params.V_D * params.V_D
          / (params.S_D * params.Tc_0 * Q2);
        d2GdPdT = 0.5 * params.V_D / (Tc * Q2);
      } else {
        G = params.S_D
          * ((temperature - Tc) * 4.0 + params.Tc_0 * 64.0 / 3.0)
          + G_disordered;
        dGdT = params.S_D * 4.0 + dGdT_disordered;
        dGdP = -params.V_D * 4.0 + dGdP_disordered;
        d2GdT2 = 0.0;
        d2GdP2 = 0.0;
        d2GdPdT = 0.0;
      }
    } else {
      Q2 = 0.0;
      G = G_disordered;
      dGdT = dGdT_disordered;
      dGdP = dGdP_disordered;
      d2GdT2 = 0.0;
      d2GdP2 = 0.0;
      d2GdPdT = 0.0;
    }
    Excesses landau_ex{
      G,
      dGdT,
      dGdP,
      d2GdT2,
      d2GdP2,
      d2GdPdT};
    // double Q = std::sqrt(Q2);
    return landau_ex;
  }

  Excesses compute_excesses(
    double pressure,
    double temperature,
    LandauHPParams params
  ) {
    double Q_0;
    double Q;
    double G;
    double dGdT;
    double dGdP;
    double d2GdT2;
    double d2GdP2;
    double d2GdPdT;
    if (params.T_0 < params.Tc_0) {
      double val = (params.Tc_0 - params.T_0) / params.Tc_0;
      Q_0 = std::sqrt(std::sqrt(val));
    } else {
      Q_0 = 0.0;
    }
    double Tc = params.Tc_0 + params.V_D
      * (pressure - params.P_0) / params.S_D;
    if (temperature < Tc) {
      double val = (Tc - temperature) / params.Tc_0;
      Q = std::sqrt(std::sqrt(val));
    } else {
      Q = 0.0;
    }
    double Q2 = Q * Q;
    double Q6 = Q2 * Q2 * Q2;
    double Q_02 = Q_0 * Q_0;
    double Q_06 = Q_02 * Q_02 * Q_02;
    G = params.Tc_0
      * params.S_D * (Q_02 - Q_06 / 3.0)
      - params.S_D * (Tc * Q2 - params.Tc_0 * Q6 / 3.0)
      - temperature * params.S_D * (Q_02 - Q2)
      + (pressure - params.P_0) * params.V_D * Q_02;
    dGdT = params.S_D * (Q2 - Q_02);
    dGdP = -params.V_D * (Q2 - Q_02);
    if (Q > constants::precision::abs_tolerance) {
      d2GdT2 = -params.S_D / (2.0 * params.Tc_0 * Q2);
      d2GdP2 = -params.V_D * params.V_D
        / (2.0 * params.S_D * params.Tc_0 * Q2);
      d2GdPdT = params.V_D / (2.0 * params.Tc_0 * Q2);
    } else {
      d2GdT2 = 0.0;
      d2GdP2 = 0.0;
      d2GdPdT = 0.0;
    }
    Excesses landau_ex{
      G,
      dGdT,
      dGdP,
      d2GdT2,
      d2GdP2,
      d2GdPdT};
    // Q;
    return landau_ex;
  }

  Excesses compute_excesses(
    double pressure,
    double temperature,
    LinearParams params
  ) {
    double G = params.delta_E
      - temperature * params.delta_S
      + pressure * params.delta_V;
    double dGdT = -params.delta_S;
    double dGdP = params.delta_V;
    double d2GdT2 = 0.0;
    double d2GdP2 = 0.0;
    double d2GdPdT = 0.0;
    Excesses linear_ex{
      G,
      dGdT,
      dGdP,
      d2GdT2,
      d2GdP2,
      d2GdPdT};
    // No Q return
    return linear_ex;
  }

  Excesses compute_excesses(
    double pressure,
    double temperature,
    BraggWilliamsParams params
  ) {
    int f_0;
    int f_1;
    if (params.factor > 0) {
      f_0 = params.factor;
      f_1 = params.factor;
    } else {
      f_0 = 1;
      f_1 = -params.factor;
    }
    double dT = 0.1;
    double dP = 1000.0;
    double G = order_gibbs(pressure, temperature, params, f_0, f_1);
    double GsubPsubT = order_gibbs(
      pressure - dP,
      temperature - dT,
      params, f_0, f_1);
    double GsubPaddT = order_gibbs(
      pressure - dP,
      temperature + dT,
      params, f_0, f_1);
    double GaddPsubT = order_gibbs(
      pressure + dP,
      temperature - dT,
      params, f_0, f_1);
    double GaddPaddT = order_gibbs(
      pressure + dP,
      temperature + dT,
      params, f_0, f_1);
    double GsubP = order_gibbs(
      pressure - dP,
      temperature,
      params, f_0, f_1);
    double GaddP = order_gibbs(
      pressure + dP,
      temperature,
      params, f_0, f_1);
    double GsubT = order_gibbs(
      pressure,
      temperature - dT,
      params, f_0, f_1);
    double GaddT = order_gibbs(
      pressure,
      temperature + dT,
      params, f_0, f_1);
    double dGdT = (GaddT - GsubT) / (2.0 * dT);
    double dGdP = (GaddP - GsubP) / (2.0 * dP);
    double d2GdT2 = (GaddT + GsubT - 2.0 * G) / (dT * dT);
    double d2GdP2 = (GaddP + GsubP - 2.0 * G) / (dP * dP);
    double d2GdPdT = (GaddPaddT - GsubPaddT - GaddPsubT + GsubPsubT)
      / (4.0 * dT * dP);
    Excesses bw_ex{
      G,
      dGdT,
      dGdP,
      d2GdT2,
      d2GdP2,
      d2GdPdT};
    // No Q return
    return bw_ex;
  }

  Excesses compute_excesses(
    double pressure,
    double temperature,
    MagneticChsParams params
  ) {
    double p = params.structural_parameter;
    double curie_T = params.curie_T_0 + pressure * params.curie_T_p;
    double curie_T2 = curie_T * curie_T;
    double tau = temperature / curie_T;
    double dtaudT = 1.0 / curie_T;
    double dtaudP = -(temperature * params.curie_T_p) / curie_T2;
    double d2taudPdT = params.curie_T_p / curie_T2;
    double d2taudP2 = 2.0 * temperature * params.curie_T_p * params.curie_T_p
      / (curie_T * curie_T2);
    double mu = params.magnetic_moment_0
      + pressure * params.magnetic_moment_p;
    double dmudP = params.magnetic_moment_p;
    double A = (518.0 / 1125.0) + (11692.0 / 15975.0) * ((1.0 / p) - 1.0);
    double f;
    double dfdtau;
    double d2fdtau2;
    if (tau < 1.0) {
      f = 1.0 - (1.0 / A)
        * (79.0 / (140.0 * p * tau)
          + (474.0 / 497.0)
          + (1.0 / p - 1.0)
          + (std::pow(tau, 3) / 6.0
            + std::pow(tau, 9) / 135.0
            + std::pow(tau, 15) / 600.0));
      dfdtau = -(1.0 / A)
        * (-79.0 / (140.0 * p * tau * tau)
          + (474.0 / 497.0)
          + (1.0 / p - 1.0)
          * (tau * tau / 2.0 + std::pow(tau, 8)
            / 15.0 + std::pow(tau, 14) / 40.0));

      d2fdtau2 = -(1.0 / A)
        * (2.0 * 79.0 / (140.0 * p * std::pow(tau, 3))
          + (474.0 / 497.0)
          * (1.0 / p - 1.0)
          * (tau
            + 8.0 * std::pow(tau, 7) / 15.0
            + 14.0 * std::pow(tau, 13) / 4.0));
    } else {
      f = -(1.0 / A)
        * (std::pow(tau, -5) / 10.0
          + std::pow(tau, -15) / 315.0
          + std::pow(tau, -25) / 1500);
      dfdtau = (1.0 / A)
        * (std::pow(tau, -6) / 2.0
          + std::pow(tau, -16) / 21.0
          + std::pow(tau, -26) / 60.0);
      d2fdtau2 = -(1.0 / A)
        * (6.0 * std::pow(tau, -7) / 2.0
          + 16.0 * std::pow(tau, -17) / 21.0
          + 26.0 * std::pow(tau, -27) / 60.0);
    }
    double dfdT = dfdtau * dtaudT;
    double d2fdT2 = d2fdtau2 * dtaudT * dtaudT;
    double dfdP = dfdtau * dtaudP;
    double d2fdP2 = d2fdtau2 * dtaudP * dtaudP + dfdtau * d2taudP2;
    double d2fdPdT = d2fdtau2 * dtaudT * dtaudP - dfdtau * d2taudPdT;
    double log_mu_p1 = std::log1p(mu);
    double G = constants::physics::gas_constant
      * temperature * log_mu_p1 * f;
    double dGdT = constants::physics::gas_constant * log_mu_p1
      * (f + temperature * dfdT);
    double d2GdT2 = constants::physics::gas_constant * log_mu_p1
      * (2.0 * dfdT + temperature * d2fdT2);
    double dGdP = constants::physics::gas_constant * temperature
      * (f * dmudP / (mu + 1.0) + dfdP * log_mu_p1);
    double val = dmudP / (mu + 1);
    double d2GdP2 = constants::physics::gas_constant * temperature
      * (-f * val * val + 2 * dfdP * val + d2fdP2 * log_mu_p1);
    double d2GdPdT = dGdP / temperature
      + (constants::physics::gas_constant
        * temperature * log_mu_p1 * d2fdPdT
        + constants::physics::gas_constant
        * temperature * dmudP / (mu + 1.0) * dfdT);
    Excesses magnetic_ex{
      G,
      dGdT,
      dGdP,
      d2GdT2,
      d2GdP2,
      d2GdPdT};
    // No Q return
    return magnetic_ex;
  }

  Excesses compute_excesses(
    double pressure [[maybe_unused]],
    double temperature,
    DebyeParams params
  ) {
    double d2GdT2;
    double f = params.Cv_inf / 3.0 / constants::physics::gas_constant;
    double G = debye::compute_helmholtz_free_energy(
      temperature, params.Theta_0, f);
    double dGdT = -debye::compute_entropy(
      temperature, params.Theta_0, f);
    double dGdP = 0.0;
    if (temperature > constants::precision::double_eps) {
      d2GdT2 = -debye::compute_molar_heat_capacity_v(
        temperature, params.Theta_0, f) / temperature;
    } else {
      d2GdT2 = 0.0;
    }
    double d2GdP2 = 0.0;
    double d2GdPdT = 0.0;
    Excesses debye_ex{
      G,
      dGdT,
      dGdP,
      d2GdT2,
      d2GdP2,
      d2GdPdT};
    // No Q return
    return debye_ex;
  }

  Excesses compute_excesses(
    double pressure [[maybe_unused]],
    double temperature,
    DebyeDeltaParams params
  ) {
    double f = params.S_inf / 3.0 / constants::physics::gas_constant;
    double G = -debye::compute_thermal_energy(
      temperature, params.Theta_0, f);
    double dGdT = -debye::compute_molar_heat_capacity_v(
      temperature, params.Theta_0, f);
    double dGdP = 0.0;
    double d2GdT2 = -debye::compute_dmolar_heat_capacity_v_dT(
      temperature, params.Theta_0, f);
    double d2GdP2 = 0.0;
    double d2GdPdT = 0.0;
    Excesses debye_ex{
      G,
      dGdT,
      dGdP,
      d2GdT2,
      d2GdP2,
      d2GdPdT};
    // No Q return
    return debye_ex;
  }

  Excesses compute_excesses(
    double pressure [[maybe_unused]],
    double temperature,
    EinsteinParams params
  ) {

    double f = params.Cv_inf / 3.0 / constants::physics::gas_constant;
    double G = einstein::compute_helmholtz_free_energy(
      temperature, params.Theta_0, f);
    double dGdT = -einstein::compute_entropy(
      temperature, params.Theta_0, f);
    double dGdP = 0.0;
    double d2GdT2;
    if (temperature > constants::precision::double_eps) {
      d2GdT2 = -einstein::compute_molar_heat_capacity_v(
        temperature, params.Theta_0, f) / temperature;
    } else {
      d2GdT2 = 0.0;
    }
    double d2GdP2 = 0.0;
    double d2GdPdT = 0.0;
    Excesses einstein_ex{
      G,
      dGdT,
      dGdP,
      d2GdT2,
      d2GdP2,
      d2GdPdT};
    // No Q return
    return einstein_ex;
  }

  Excesses compute_excesses(
    double pressure [[maybe_unused]],
    double temperature,
    EinsteinDeltaParams params
  ) {
    double f = params.S_inf / 3.0 / constants::physics::gas_constant;
    double G = -einstein::compute_thermal_energy(
      temperature, params.Theta_0, f);
    double dGdT = -einstein::compute_molar_heat_capacity_v(
      temperature, params.Theta_0, f);
    double dGdP = 0.0;
    double d2GdT2 = -einstein::compute_dmolar_heat_capacity_v_dT(
      temperature, params.Theta_0, f);
    double d2GdP2 = 0.0;
    double d2GdPdT = 0.0;
    Excesses einstein_ex{
      G,
      dGdT,
      dGdP,
      d2GdT2,
      d2GdP2,
      d2GdPdT};
    // No Q return
    return einstein_ex;
  }

} // End namespace excesses
