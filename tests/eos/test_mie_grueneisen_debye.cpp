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
#include <cmath>

using namespace Catch::Matchers;

TEST_CASE("Test validate parameters", "[mgd][eos]") {
  MGD3 mgd;
  SECTION("Missing V_0") {
    MineralParams params;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing K_0") {
    MineralParams params;
    params.V_0 = 11.24e-6;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing Kprime_0") {
    MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing molar mass") {
    MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing napfu") {
    MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing debye_0") {
    MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing grueneisen_0") {
    MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing q_0") {
    MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }

  SECTION("Optional Parameters") {
    MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_FALSE(params.T_0.has_value());
    REQUIRE_FALSE(params.P_0.has_value());
    REQUIRE_FALSE(params.E_0.has_value());
    REQUIRE_FALSE(params.F_0.has_value());
    REQUIRE_FALSE(params.G_0.has_value());
    REQUIRE_FALSE(params.Gprime_0.has_value());
    mgd.validate_parameters(params);
    REQUIRE(*params.T_0 == 300.0);
    REQUIRE(*params.P_0 == 0);
    REQUIRE(*params.E_0 == 0);
    REQUIRE(*params.F_0 == 0);
    REQUIRE(std::isnan(*params.G_0));
    REQUIRE(std::isnan(*params.Gprime_0));
  }
  // TODO: Check warnings
}

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
    //CHECK_THAT(mgd3.compute_debye_temperature(1.0, params),
    //  WithinRel(*params.debye_0, tol_rel) ||
    //  WithinAbs(*params.debye_0, tol_abs));
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

TEST_CASE("Shear modulus expansion", "[birch-murnaghan][eos]") {
  // Set up test P, T, V
  double P = 25.e9;
  double T = 2000.0;
  double V = 10.0e-6;
  // Set up test params
  MineralParams params;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.debye_0 = 773.0;
  params.grueneisen_0 = 1.5;
  params.q_0 = 1.5;
  params.molar_mass = 0.0403;
  params.napfu = 2;
  MGD2 mgd2;
  MGD3 mgd3;
  // Validating should set G0/G' to nan
  mgd3.validate_parameters(params);
  REQUIRE_NOTHROW(mgd2.compute_shear_modulus(P, T, V, params));
  REQUIRE_NOTHROW(mgd3.compute_shear_modulus(P, T, V, params));
  CHECK(std::isnan(mgd2.compute_shear_modulus(P, T, V, params)));
  CHECK(std::isnan(mgd3.compute_shear_modulus(P, T, V, params)));
  // Set G0/G'0
  params.G_0 = 131.0e9;
  params.Gprime_0 = 2.1;
  CHECK_FALSE(std::isnan(mgd3.compute_shear_modulus(P, T, V, params)));
  // Ensure expansion order different
  CHECK_FALSE(
    mgd2.compute_shear_modulus(P, T, V, params) ==
    mgd3.compute_shear_modulus(P, T, V, params)
  );
}