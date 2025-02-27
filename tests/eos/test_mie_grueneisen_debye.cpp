/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "tolerances.hpp"
#include "burnman/eos/mie_grueneisen_debye.hpp"

using namespace Catch::Matchers;

TEST_CASE("Check reference conditions", "[mgd][eos]") {
  // Set up test params
  MineralParams params;
  params.E_0 = 0.0;
  params.T_0 = 300.0;
  params.P_0 = 0.0;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.G_0 = 131.0e9;
  params.Gprime_0 = 2.1;
  params.debye_0 = 773.0;
  params.grueneisen_0 = 1.5;
  params.q_0 = 1.5;
  params.molar_mass = 0.0403;
  params.napfu = 2;

  MGD2 mgd2;
  MGD3 mgd3;
  REQUIRE_NOTHROW(mgd3.validate_parameters(params));

  SECTION("Independent functions") {
    double T = *params.T_0;
    double P = *params.P_0;
    double V = *params.V_0;
    CHECK_THAT(mgd3.compute_pressure(T, V, params),
      WithinRel(*params.P_0, tol_rel) || WithinAbs(*params.P_0, tol_abs));
    //CHECK_THAT(mgd3.compute_volume(P, T, params),
    //  WithinRel(*params.V_0, tol_rel) || WithinAbs(*params.V_0, tol_abs));
  }
  SECTION("V dependent functions") {
    auto P = GENERATE(0.0, 10.0, 25.e9);
    auto T = GENERATE(300.0, 2000.0);
    double V = *params.V_0;    
    CHECK_THAT(mgd3.compute_grueneisen_parameter(P, T, V, params),
      WithinRel(*params.grueneisen_0, tol_rel) ||
      WithinAbs(*params.grueneisen_0, tol_abs));
    // x = V_0 / V
    CHECK_THAT(mgd3.compute_debye_temperature(1.0, params),
      WithinRel(*params.debye_0, tol_rel) ||
      WithinAbs(*params.debye_0, tol_abs));
  }
  SECTION("V-T dependent functions") {
    auto P = GENERATE(0.0, 10.0, 25.e9);
    double T = *params.T_0;
    double V = *params.V_0;
    CHECK_THAT(mgd3.compute_isothermal_bulk_modulus_reuss(P, T, V, params),
      WithinRel(*params.K_0, tol_rel) || WithinAbs(*params.K_0, tol_abs));
    CHECK_THAT(mgd3.compute_shear_modulus(P, T, V, params),
      WithinRel(*params.G_0, tol_rel) || WithinAbs(*params.G_0, tol_abs));
    CHECK_THAT(mgd2.compute_shear_modulus(P, T, V, params),
      WithinRel(*params.G_0, tol_rel) || WithinAbs(*params.G_0, tol_abs));
    CHECK_THAT(mgd3.compute_helmholtz_free_energy(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
  }
  SECTION("P-T-V dependent functions") {
    MineralParams params_b = params;
    params_b.T_0 = 500;
    params_b.P_0 = 1.e6;
    double P = *params.P_0;
    double T = *params.T_0;
    double V = *params.V_0;
    double P2 = *params_b.P_0;
    double T2 = *params_b.T_0;
    CHECK_THAT(mgd3.compute_gibbs_free_energy(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
    CHECK_THAT(mgd3.compute_gibbs_free_energy(P2, T2, V, params_b),
      WithinRel(P2*V, tol_rel) || WithinAbs(P2*V, tol_abs));
  }
  SECTION("Zero T known functions") {
    auto P = GENERATE(0.0, 10.0, 25.e9);
    double T = 0.0;
    double V = *params.V_0;
    CHECK_THAT(mgd3.compute_molar_heat_capacity_v(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
    CHECK_THAT(mgd3.compute_molar_heat_capacity_p(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
    CHECK_THAT(mgd3.compute_thermal_expansivity(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
    // At T_ref: K_T = K_0; At T=0: K_S = K_T
    CHECK_THAT(mgd3.compute_isentropic_bulk_modulus_reuss(P, T, V, params),
      !WithinRel(*params.K_0, tol_rel) && !WithinAbs(*params.K_0, tol_abs));
    // Copy params to set T_0 = 0
    MineralParams params_T0 = params;
    params_T0.T_0 = 0;
    CHECK_THAT(mgd3.compute_isentropic_bulk_modulus_reuss(P, T, V, params_T0),
      WithinRel(*params.K_0, tol_rel) || WithinAbs(*params.K_0, tol_abs));
  }
}
  
