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
#include <tuple>
#include <map>

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

TEST_CASE("Shear modulus expansion", "[mgd][eos]") {
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


TEST_CASE("MGD python reference values", "[mgd][eos]") {
  // Set up test params
  MineralParams params;
  params.T_0 = 300.0;
  params.P_0 = 0.0;
  params.E_0 = 0.0;
  params.F_0 = 0.0;
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

  SECTION("Test volume") {
    MGD3 mgd;
    double T_a = 800.0;
    double T_b = 2000.0;
    double V_a = 0.9 * (*params.V_0);
    double V_b = 0.5 * (*params.V_0);
    double P_aa = mgd.compute_pressure(T_a, V_a, params);
    double P_ab = mgd.compute_pressure(T_a, V_b, params);
    double P_ba = mgd.compute_pressure(T_b, V_a, params);
    double P_bb = mgd.compute_pressure(T_b, V_b, params);
    CHECK_THAT(mgd.compute_volume(P_aa, T_a, params),
      WithinRel(V_a, tol_rel) || WithinAbs(V_a, tol_abs));
    CHECK_THAT(mgd.compute_volume(P_ab, T_a, params),
      WithinRel(V_b, tol_rel) || WithinAbs(V_b, tol_abs));
    CHECK_THAT(mgd.compute_volume(P_ba, T_b, params),
      WithinRel(V_a, tol_rel) || WithinAbs(V_a, tol_abs));
    CHECK_THAT(mgd.compute_volume(P_bb, T_b, params),
      WithinRel(V_b, tol_rel) || WithinAbs(V_b, tol_abs));
  }
  SECTION("V dependent functions") {
    struct TestData {
      double input;
      double expected_gamma;
    };
    MGD3 mgd3;
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
    CHECK_THAT(bm3.compute_grueneisen_parameter(
        P, T, test_data.input*(*params.V_0), params
      ),
      WithinRel(test_data.expected_gamma, tol_rel) ||
      WithinAbs(test_data.expected_gamma, tol_abs));
  }
  SECTION("T-V dependent functions") {
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
    // Generate combos
    auto T = GENERATE(T1, T2, T3);
    auto x = GENERATE(x1, x2, x3);
    CAPTURE[T];
    CAPTURE[x];
    double V = *params.V_0 * x;
    auto key = std::make_tuple(T, x);
    //
    CHECK_THAT(mgd3.compute_pressure(
        T, V, params
      ),
      WithinRel(ref_P[key], tol_rel) ||
      WithinAbs(ref_P[key], tol_abs));
    CHECK_THAT(mgd3.compute_isothermal_bulk_modulus_reuss(
        P, T, V, params
      ),
      WithinRel(ref_KT[key], tol_rel) ||
      WithinAbs(ref_KT[key], tol_abs));
    CHECK_THAT(mgd3.compute_isentropic_bulk_modulus_reuss(
        P, T, V, params
      ),
      WithinRel(ref_KS[key], tol_rel) ||
      WithinAbs(ref_KS[key], tol_abs));
    CHECK_THAT(mgd2.compute_shear_modulus(
        P, T, V, params
      ),
      WithinRel(ref_G2[key], tol_rel) ||
      WithinAbs(ref_G2[key], tol_abs));
    CHECK_THAT(mgd3.compute_shear_modulus(
        P, T, V, params
      ),
      WithinRel(ref_G3[key], tol_rel) ||
      WithinAbs(ref_G3[key], tol_abs));
    CHECK_THAT(mgd3.compute_molar_heat_capacity_v(
        P, T, V, params
      ),
      WithinRel(ref_Cv[key], tol_rel) ||
      WithinAbs(ref_Cv[key], tol_abs));
    CHECK_THAT(mgd3.compute_molar_heat_capacity_p(
        P, T, V, params
      ),
      WithinRel(ref_Cp[key], tol_rel) ||
      WithinAbs(ref_Cp[key], tol_abs));
    CHECK_THAT(mgd3.compute_entropy(
        P, T, V, params
      ),
      WithinRel(ref_S[key], tol_rel) ||
      WithinAbs(ref_S[key], tol_abs));
    CHECK_THAT(mgd3.compute_helmholtz_free_energy(
        P, T, V, params
      ),
      WithinRel(ref_F[key], tol_rel) ||
      WithinAbs(ref_F[key], tol_abs));
    CHECK_THAT(mgd3.compute_molar_internal_energy(
        P, T, V, params
      ),
      WithinRel(ref_E[key], tol_rel) ||
      WithinAbs(ref_E[key], tol_abs));
    CHECK_THAT(mgd3.compute_thermal_expansivity(
        P, T, V, params
      ),
      WithinRel(ref_alpha[key], tol_rel) ||
      WithinAbs(ref_alpha[key], tol_abs));
  }
  SECTION("P-T-V dependent functions") {
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
    auto P = GENERATE(P1, P2);
    auto T = GENERATE(T1, T2);
    auto x = GENERATE(x1, x2);
    CAPTURE[P];
    CAPTURE[T];
    CAPTURE[x];
    double V = *params.V_0 * x;
    auto key = std::make_tuple(P, T, x);
    //
    CHECK_THAT(mgd3.compute_gibbs_free_energy(
        P, T, V, params
      ),
      WithinRel(ref_G[key], tol_rel) ||
      WithinAbs(ref_G[key], tol_abs));
    CHECK_THAT(mgd3.compute_enthalpy(
        P, T, V, params
      ),
      WithinRel(ref_H[key], tol_rel) ||
      WithinAbs(ref_H[key], tol_abs));
  }
}
