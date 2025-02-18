/*
  TODO: Copyright Notice!
*/
#include <cmath>
#include "burnman/eos/property_modifiers.hpp"
#include "burnman/eos/debye.hpp"
#include "burnman/eos/einstein.hpp"
#include "burnman/utils/constants.hpp"

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
      double Q2 = std::sqrt((Tc - temperature) / params.Tc_0);
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
    G = params.Tc_0 * params.S_D
      * (Q_02 - Q_06 / 3.0)
      - temperature
      * params.S_D * (Q_02 - Q2)
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
    BraggWilliamsParams params);

  Excesses compute_excesses(
    double pressure,
    double temperature,
    MagneticChsParams params);

  Excesses compute_excesses(
    double pressure,
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
    double pressure,
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
    double pressure,
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
    double pressure,
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
