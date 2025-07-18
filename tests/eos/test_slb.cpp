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
#include "burnman/utils/eos.hpp"
#include "burnman/eos/slb.hpp"
#include <cmath>
#include <tuple>
#include <map>
using namespace Catch::Matchers;

TEST_CASE("Test validate parameters", "[eos][slb]") {
  SLB3 slb;
  SECTION("Missing V_0") {
    MineralParams params;
    params.molar_mass = 0.0055;
    params.napfu = 1;
    params.K_0 = 1.6e11;
    params.Kprime_0 = 6.0;
    params.F_0 = 11.82;
    params.debye_0 = 398.0;
    params.grueneisen_0 = 1.7;
    params.q_0 = 0.9;
    params.G_0 = 8.1e10;
    params.Gprime_0 = 1.9;
    params.eta_s_0 = 5.9;
    params.bel_0 = 0.004;
    params.gel = 1.5;
    REQUIRE_THROWS(slb.validate_parameters(params));
  }
  SECTION("Missing K_0") {
    MineralParams params;
    params.molar_mass = 0.0055;
    params.napfu = 1;
    params.V_0 = 7.0e-06;
    params.Kprime_0 = 6.0;
    params.F_0 = 11.82;
    params.debye_0 = 398.0;
    params.grueneisen_0 = 1.7;
    params.q_0 = 0.9;
    params.G_0 = 8.1e10;
    params.Gprime_0 = 1.9;
    params.eta_s_0 = 5.9;
    params.bel_0 = 0.004;
    params.gel = 1.5;
    REQUIRE_THROWS(slb.validate_parameters(params));
  }
  SECTION("Missing Kprime_0") {
    MineralParams params;
    params.molar_mass = 0.0055;
    params.napfu = 1;
    params.V_0 = 7.0e-06;
    params.K_0 = 1.6e11;
    params.F_0 = 11.82;
    params.debye_0 = 398.0;
    params.grueneisen_0 = 1.7;
    params.q_0 = 0.9;
    params.G_0 = 8.1e10;
    params.Gprime_0 = 1.9;
    params.eta_s_0 = 5.9;
    params.bel_0 = 0.004;
    params.gel = 1.5;
    REQUIRE_THROWS(slb.validate_parameters(params));
  }
  SECTION("Missing molar mass") {
    MineralParams params;
    params.napfu = 1;
    params.V_0 = 7.0e-06;
    params.K_0 = 1.6e11;
    params.Kprime_0 = 6.0;
    params.F_0 = 11.82;
    params.debye_0 = 398.0;
    params.grueneisen_0 = 1.7;
    params.q_0 = 0.9;
    params.G_0 = 8.1e10;
    params.Gprime_0 = 1.9;
    params.eta_s_0 = 5.9;
    params.bel_0 = 0.004;
    params.gel = 1.5;
    REQUIRE_THROWS(slb.validate_parameters(params));
  }
  SECTION("Missing napfu") {
    MineralParams params;
    params.molar_mass = 0.0055;
    params.V_0 = 7.0e-06;
    params.K_0 = 1.6e11;
    params.Kprime_0 = 6.0;
    params.F_0 = 11.82;
    params.debye_0 = 398.0;
    params.grueneisen_0 = 1.7;
    params.q_0 = 0.9;
    params.G_0 = 8.1e10;
    params.Gprime_0 = 1.9;
    params.eta_s_0 = 5.9;
    params.bel_0 = 0.004;
    params.gel = 1.5;
    REQUIRE_THROWS(slb.validate_parameters(params));
  }
  SECTION("Missing debye_0") {
    MineralParams params;
    params.molar_mass = 0.0055;
    params.napfu = 1;
    params.V_0 = 7.0e-06;
    params.K_0 = 1.6e11;
    params.Kprime_0 = 6.0;
    params.F_0 = 11.82;
    params.grueneisen_0 = 1.7;
    params.q_0 = 0.9;
    params.G_0 = 8.1e10;
    params.Gprime_0 = 1.9;
    params.eta_s_0 = 5.9;
    params.bel_0 = 0.004;
    params.gel = 1.5;
    REQUIRE_THROWS(slb.validate_parameters(params));
  }
  SECTION("Missing grueneisen_0") {
    MineralParams params;
    params.molar_mass = 0.0055;
    params.napfu = 1;
    params.V_0 = 7.0e-06;
    params.K_0 = 1.6e11;
    params.Kprime_0 = 6.0;
    params.F_0 = 11.82;
    params.debye_0 = 398.0;
    params.q_0 = 0.9;
    params.G_0 = 8.1e10;
    params.Gprime_0 = 1.9;
    params.eta_s_0 = 5.9;
    params.bel_0 = 0.004;
    params.gel = 1.5;
    REQUIRE_THROWS(slb.validate_parameters(params));
  }
  SECTION("Missing q_0") {
    MineralParams params;
    params.molar_mass = 0.0055;
    params.napfu = 1;
    params.V_0 = 7.0e-06;
    params.K_0 = 1.6e11;
    params.Kprime_0 = 6.0;
    params.F_0 = 11.82;
    params.debye_0 = 398.0;
    params.grueneisen_0 = 1.7;
    params.G_0 = 8.1e10;
    params.Gprime_0 = 1.9;
    params.eta_s_0 = 5.9;
    params.bel_0 = 0.004;
    params.gel = 1.5;
    REQUIRE_THROWS(slb.validate_parameters(params));
  }
  SECTION("SLB3Conductive; Missing bel_0") {
    SLB3Conductive slb_c;
    MineralParams params;
    params.molar_mass = 0.0055;
    params.napfu = 1;
    params.V_0 = 7.0e-06;
    params.K_0 = 1.6e11;
    params.Kprime_0 = 6.0;
    params.F_0 = 11.82;
    params.debye_0 = 398.0;
    params.grueneisen_0 = 1.7;
    params.q_0 = 0.9;
    params.G_0 = 8.1e10;
    params.Gprime_0 = 1.9;
    params.eta_s_0 = 5.9;
    params.gel = 1.5;
    REQUIRE_THROWS(slb_c.validate_parameters(params));
  }
  SECTION("SLB3Conductive; Missing gel") {
    SLB3Conductive slb_c;
    MineralParams params;
    params.molar_mass = 0.0055;
    params.napfu = 1;
    params.V_0 = 7.0e-06;
    params.K_0 = 1.6e11;
    params.Kprime_0 = 6.0;
    params.F_0 = 11.82;
    params.debye_0 = 398.0;
    params.grueneisen_0 = 1.7;
    params.q_0 = 0.9;
    params.G_0 = 8.1e10;
    params.Gprime_0 = 1.9;
    params.eta_s_0 = 5.9;
    params.bel_0 = 0.004;
    REQUIRE_THROWS(slb_c.validate_parameters(params));
  }

  SECTION("Optional Parameters") {
    MineralParams params;
    params.molar_mass = 0.0055;
    params.napfu = 1;
    params.V_0 = 7.0e-06;
    params.K_0 = 1.6e11;
    params.Kprime_0 = 6.0;
    params.debye_0 = 398.0;
    params.grueneisen_0 = 1.7;
    params.q_0 = 0.9;
    REQUIRE_FALSE(params.T_0.has_value());
    REQUIRE_FALSE(params.P_0.has_value());
    REQUIRE_FALSE(params.E_0.has_value());
    REQUIRE_FALSE(params.F_0.has_value());
    REQUIRE_FALSE(params.G_0.has_value());
    REQUIRE_FALSE(params.Gprime_0.has_value());
    REQUIRE_FALSE(params.eta_s_0.has_value());
    REQUIRE_NOTHROW(slb.validate_parameters(params));
    REQUIRE(*params.T_0 == 300.0);
    REQUIRE(*params.P_0 == 0);
    REQUIRE(*params.E_0 == 0);
    REQUIRE(*params.F_0 == 0);
    REQUIRE(std::isnan(*params.G_0));
    REQUIRE(std::isnan(*params.Gprime_0));
    REQUIRE(std::isnan(*params.eta_s_0));
  }
}

TEST_CASE("Check reference conditions", "[eos][slb]") {
  // Set up test params
  MineralParams params;
  params.E_0 = 0.0;
  params.F_0 = -2.05e6;
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
  params.eta_s_0 = 4.5;
  params.molar_mass = 0.04;
  params.napfu = 2;
  SLB2 slb2;
  SLB3 slb3;
  SLB3Conductive slb3_c;
  REQUIRE_NOTHROW(slb3.validate_parameters(params));

  SECTION("P-T-V dependent functions") {
    double T = *params.T_0;
    double P = *params.P_0;
    double V = *params.V_0;
    CHECK_THAT(slb3.compute_pressure(T, V, params),
      WithinRel(*params.P_0, tol_rel) || WithinAbs(*params.P_0, tol_abs));
    CHECK_THAT(slb3.compute_volume(P, T, params),
      WithinRel(*params.V_0, tol_rel) || WithinAbs(*params.V_0, tol_abs));
    CHECK_THAT(slb3_c.compute_pressure(T, V, params),
      WithinRel(*params.P_0, tol_rel) || WithinAbs(*params.P_0, tol_abs));
    CHECK_THAT(slb3_c.compute_volume(P, T, params),
      WithinRel(*params.V_0, tol_rel) || WithinAbs(*params.V_0, tol_abs));
    CHECK_THAT(slb3.compute_gibbs_free_energy(P, T, V, params),
      WithinRel(*params.F_0, tol_rel) || WithinAbs(*params.F_0, tol_abs));
  }
  SECTION("V dependent functions") {
    auto P = GENERATE(0.0, 10.0, 25.e9);
    auto T = GENERATE(300.0, 2000.0);
    double V = *params.V_0;
    CHECK_THAT(slb3.compute_grueneisen_parameter(P, T, V, params),
      WithinRel(*params.grueneisen_0, tol_rel) ||
      WithinAbs(*params.grueneisen_0, tol_abs));
  }
  SECTION("V-T dependent functions") {
    auto P = GENERATE(0.0, 10.0, 25.e9);
    double T = *params.T_0;
    double V = *params.V_0;
    CHECK_THAT(slb3.compute_isothermal_bulk_modulus_reuss(P, T, V, params),
      WithinRel(*params.K_0, tol_rel) || WithinAbs(*params.K_0, tol_abs));
    CHECK_THAT(slb3.compute_shear_modulus(P, T, V, params),
      WithinRel(*params.G_0, tol_rel) || WithinAbs(*params.G_0, tol_abs));
    CHECK_THAT(slb3.compute_shear_modulus(P, T, V, params),
      WithinRel(*params.G_0, tol_rel) || WithinAbs(*params.G_0, tol_abs));
    CHECK_THAT(slb3.compute_helmholtz_free_energy(P, T, V, params),
      WithinRel(*params.F_0, tol_rel) || WithinAbs(*params.F_0, tol_abs));
  }
  SECTION("Zero T known functions") {
    auto P = GENERATE(0.0, 10.0, 25.e9);
    double T = 0.0;
    double V = *params.V_0;
    CHECK_THAT(slb3.compute_molar_heat_capacity_v(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
    CHECK_THAT(slb3.compute_molar_heat_capacity_p(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
    CHECK_THAT(slb3.compute_thermal_expansivity(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
    // At T_ref: K_T = K_0; At T=0: K_S = K_T
    CHECK_THAT(slb3.compute_isentropic_bulk_modulus_reuss(P, T, V, params),
      !WithinRel(*params.K_0, tol_rel) && !WithinAbs(*params.K_0, tol_abs));
    // Copy params to set T_0 = 0
    MineralParams params_T0 = params;
    params_T0.T_0 = 0;
    CHECK_THAT(slb3.compute_isentropic_bulk_modulus_reuss(P, T, V, params_T0),
      WithinRel(*params.K_0, tol_rel) || WithinAbs(*params.K_0, tol_abs));
  }
}

TEST_CASE("Shear modulus expansion", "[eos][slb]") {
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
  SLB2 slb2;
  SLB3 slb3;
  // Validating should set G0/G' to nan
  slb3.validate_parameters(params);
  REQUIRE_NOTHROW(slb2.compute_shear_modulus(P, T, V, params));
  REQUIRE_NOTHROW(slb3.compute_shear_modulus(P, T, V, params));
  CHECK(std::isnan(slb2.compute_shear_modulus(P, T, V, params)));
  CHECK(std::isnan(slb3.compute_shear_modulus(P, T, V, params)));
  // Set G0/G'0
  params.G_0 = 131.0e9;
  params.Gprime_0 = 2.1;
  CHECK_FALSE(std::isnan(slb2.compute_shear_modulus(P, T, V, params)));
  CHECK_FALSE(std::isnan(slb3.compute_shear_modulus(P, T, V, params)));
  // Ensure expansion order different
  double G_slb2 = slb2.compute_shear_modulus(P, T, V, params);
  double G_slb3 = slb3.compute_shear_modulus(P, T, V, params);
  CHECK_THAT(G_slb2,
    !WithinRel(G_slb3, 1e-5) && !WithinAbs(G_slb3, 1e-5));
}

TEST_CASE("Test volume", "[eos][slb]") {
  // Set up test params
  MineralParams params;
  params.T_0 = 300.0;
  params.P_0 = 0.0;
  params.E_0 = 0.0;
  params.F_0 = -1e5;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.G_0 = 131.0e9;
  params.Gprime_0 = 2.1;
  params.debye_0 = 773.0;
  params.grueneisen_0 = 1.5;
  params.q_0 = 1.4;
  params.molar_mass = 0.0403;
  params.napfu = 2;
  SLB3 slb;
  double T_a = 800.0;
  double T_b = 2000.0;
  double V_a = 0.9 * (*params.V_0);
  double V_b = 0.5 * (*params.V_0);
  double P_aa = slb.compute_pressure(T_a, V_a, params);
  double P_ab = slb.compute_pressure(T_a, V_b, params);
  double P_ba = slb.compute_pressure(T_b, V_a, params);
  double P_bb = slb.compute_pressure(T_b, V_b, params);
  CHECK_THAT(slb.compute_volume(P_aa, T_a, params),
    WithinRel(V_a, tol_rel) || WithinAbs(V_a, tol_abs));
  CHECK_THAT(slb.compute_volume(P_ab, T_a, params),
    WithinRel(V_b, tol_rel) || WithinAbs(V_b, tol_abs));
  CHECK_THAT(slb.compute_volume(P_ba, T_b, params),
    WithinRel(V_a, tol_rel) || WithinAbs(V_a, tol_abs));
  CHECK_THAT(slb.compute_volume(P_bb, T_b, params),
    WithinRel(V_b, tol_rel) || WithinAbs(V_b, tol_abs));
}

TEST_CASE("SLB python reference values", "[eos][slb]") {
  // Set up test params
  MineralParams params;
  params.T_0 = 300.0;
  params.P_0 = 0.0;
  params.E_0 = 0.0;
  params.F_0 = -1368283.0;
  params.V_0 = 2.4445e-5;
  params.K_0 = 2.505e11;
  params.Kprime_0 = 4.14;
  params.G_0 = 1.729e11;
  params.Gprime_0 = 1.69;
  params.debye_0 = 734.0;
  params.grueneisen_0 = 1.565;
  params.q_0 = 1.10945;
  params.eta_s_0 = 2.565;
  params.bel_0 = 0.00411;
  params.gel = 1.476;
  params.molar_mass = 0.1003887;
  params.napfu = 5;

  SECTION("V dependent functions") {
    struct TestData {
      double input;
      double expected_gamma;
    };
    SLB3 slb3;
    // P & T unused
    double P = 2.0e9;
    double T = 500.0;
    // TODO: Data
    auto test_data = GENERATE(
      TestData{0.99, },
      TestData{0.98, },
      TestData{0.95, },
      TestData{0.80, },
      TestData{0.40, }
    );
    CAPTURE(test_data.input);
    CHECK_THAT(slb3.compute_grueneisen_parameter(
        P, T, test_data.input*(*params.V_0), params
      ),
      WithinRel(test_data.expected_gamma, tol_rel) ||
      WithinAbs(test_data.expected_gamma, tol_abs));
  }
  SECTION("T-V dependent functions") {
    SLB3 slb3;
    SLB2 slb2;
    SLB3Conductive slb3_c;
    double P = 2.e9;
    double T1 = 300, T2 = 800, T3 = 2500;
    double x1 = 0.99, x2 = 0.80, x3 = 0.40;
    // Reference data
    std::map<std::tuple<double, double>, double> ref_P = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_P_c = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_KT = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_KT_c = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_KS = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_G2 = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_G3 = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_Cv = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_Cv_c = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_Cp = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_S = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_S_c = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_F = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_F_c = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_E = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_alpha = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_alpha_c = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    std::map<std::tuple<double, double>, double> ref_gamma_c = {
      {{T1, x1}, },
      {{T1, x2}, },
      {{T1, x3}, },
      {{T2, x1}, },
      {{T2, x2}, },
      {{T2, x3}, },
      {{T3, x1}, },
      {{T3, x2}, },
      {{T3, x3}, }
    };
    // Generate combos
    // double T1 = 300, T2 = 800, T3 = 2500;
    // double x1 = 0.99, x2 = 0.80, x3 = 0.40;
    auto T = GENERATE(300.0, 800.0, 2500.0);
    auto x = GENERATE(0.99, 0.8, 0.4);
    CAPTURE(T);
    CAPTURE(x);
    double V = *params.V_0 * x;
    auto key = std::make_tuple(T, x);
    CHECK_THAT(slb3.compute_pressure(
        T, V, params
      ),
      WithinRel(ref_P[key], tol_rel) ||
      WithinAbs(ref_P[key], tol_abs));
    CHECK_THAT(slb3_c.compute_pressure(
        T, V, params
      ),
      WithinRel(ref_P_c[key], tol_rel) ||
      WithinAbs(ref_P_c[key], tol_abs));
    CHECK_THAT(slb3.compute_isothermal_bulk_modulus_reuss(
        P, T, V, params
      ),
      WithinRel(ref_KT[key], tol_rel) ||
      WithinAbs(ref_KT[key], tol_abs));
    CHECK_THAT(slb3_c.compute_isothermal_bulk_modulus_reuss(
        P, T, V, params
      ),
      WithinRel(ref_KT_c[key], tol_rel) ||
      WithinAbs(ref_KT_c[key], tol_abs));
    CHECK_THAT(slb3.compute_isentropic_bulk_modulus_reuss(
        P, T, V, params
      ),
      WithinRel(ref_KS[key], tol_rel) ||
      WithinAbs(ref_KS[key], tol_abs));
    CHECK_THAT(slb2.compute_shear_modulus(
        P, T, V, params
      ),
      WithinRel(ref_G2[key], tol_rel) ||
      WithinAbs(ref_G2[key], tol_abs));
    CHECK_THAT(slb3.compute_shear_modulus(
        P, T, V, params
      ),
      WithinRel(ref_G3[key], tol_rel) ||
      WithinAbs(ref_G3[key], tol_abs));
    CHECK_THAT(slb3.compute_molar_heat_capacity_v(
        P, T, V, params
      ),
      WithinRel(ref_Cv[key], tol_rel) ||
      WithinAbs(ref_Cv[key], tol_abs));
    CHECK_THAT(slb3_c.compute_molar_heat_capacity_v(
        P, T, V, params
      ),
      WithinRel(ref_Cv_c[key], tol_rel) ||
      WithinAbs(ref_Cv_c[key], tol_abs));
    CHECK_THAT(slb3.compute_molar_heat_capacity_p(
        P, T, V, params
      ),
      WithinRel(ref_Cp[key], tol_rel) ||
      WithinAbs(ref_Cp[key], tol_abs));
    CHECK_THAT(slb3.compute_entropy(
        P, T, V, params
      ),
      WithinRel(ref_S[key], tol_rel) ||
      WithinAbs(ref_S[key], tol_abs));
    CHECK_THAT(slb3_c.compute_entropy(
        P, T, V, params
      ),
      WithinRel(ref_S_c[key], tol_rel) ||
      WithinAbs(ref_S_c[key], tol_abs));
    CHECK_THAT(slb3.compute_helmholtz_free_energy(
        P, T, V, params
      ),
      WithinRel(ref_F[key], tol_rel) ||
      WithinAbs(ref_F[key], tol_abs));
    CHECK_THAT(slb3_c.compute_helmholtz_free_energy(
        P, T, V, params
      ),
      WithinRel(ref_F_c[key], tol_rel) ||
      WithinAbs(ref_F_c[key], tol_abs));
    CHECK_THAT(slb3.compute_molar_internal_energy(
        P, T, V, params
      ),
      WithinRel(ref_E[key], tol_rel) ||
      WithinAbs(ref_E[key], tol_abs));

    CHECK_THAT(slb3.compute_thermal_expansivity(
        P, T, V, params
      ),
      WithinRel(ref_alpha[key], tol_rel) ||
      WithinAbs(ref_alpha[key], tol_abs));
    CHECK_THAT(slb3_c.compute_thermal_expansivity(
        P, T, V, params
      ),
      WithinRel(ref_alpha_c[key], tol_rel) ||
      WithinAbs(ref_alpha_c[key], tol_abs));
    CHECK_THAT(slb3_c.compute_grueneisen_parameter(
        P, T, V, params
      ),
      WithinRel(ref_gamma_c[key], tol_rel) ||
      WithinAbs(ref_gamma_c[key], tol_abs));
  }
  SECTION("P-T-V dependent functions") {
    SLB3 slb3;
    double P1 = 1.e9, P2 = 54.e9;
    double T1 = 800, T2 = 2500;
    double x1 = 0.99, x2 = 0.4;
    std::map<std::tuple<double, double, double>, double> ref_G = {
        {{P1, T1, x1}, },
        {{P1, T1, x2}, },
        {{P1, T2, x1}, },
        {{P1, T2, x2}, },
        {{P2, T1, x1}, },
        {{P2, T1, x2}, },
        {{P2, T2, x1}, },
        {{P2, T2, x2}, }
    };
    std::map<std::tuple<double, double, double>, double> ref_H = {
        {{P1, T1, x1}, },
        {{P1, T1, x2}, },
        {{P1, T2, x1}, },
        {{P1, T2, x2}, },
        {{P2, T1, x1}, },
        {{P2, T1, x2}, },
        {{P2, T2, x1}, },
        {{P2, T2, x2}, }
    };
    // Generate combos
    auto P = GENERATE(1.e9, 54.e9);
    auto T = GENERATE(800.0, 2500.0);
    auto x = GENERATE(0.99, 0.4);
    CAPTURE(P);
    CAPTURE(T);
    CAPTURE(x);
    double V = *params.V_0 * x;
    auto key = std::make_tuple(P, T, x);
    CHECK_THAT(slb3.compute_gibbs_free_energy(
        P, T, V, params
      ),
      WithinRel(ref_G[key], tol_rel) ||
      WithinAbs(ref_G[key], tol_abs));
    CHECK_THAT(slb3.compute_enthalpy(
        P, T, V, params
      ),
      WithinRel(ref_H[key], tol_rel) ||
      WithinAbs(ref_H[key], tol_abs));
  }
}
