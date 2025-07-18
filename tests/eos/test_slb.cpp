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
      TestData{0.99, 1.547842795607494},
      TestData{0.98, 1.5310862961658358},
      TestData{0.95, 1.4830708730153435},
      TestData{0.80, 1.2835471880685585},
      TestData{0.40, 0.9374733906527063}
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
      {{T1, x1}, 2570544830.23939},
      {{T1, x2}, 88892646063.63199},
      {{T1, x3}, 1585790219388.1765},
      {{T2, x1}, 6143041512.372141},
      {{T2, x2}, 92267463770.23813},
      {{T2, x3}, 1588778011889.625},
      {{T3, x1}, 19517827295.104595},
      {{T3, x2}, 105840760451.95923},
      {{T3, x3}, 1607063207497.6838}
    };
    std::map<std::tuple<double, double>, double> ref_P_c = {
      {{T1, x1}, 2570544830.23939},
      {{T1, x2}, 88892646063.63199},
      {{T1, x3}, 1585790219388.1765},
      {{T2, x1}, 6210960809.706225},
      {{T2, x2}, 92328831727.23389},
      {{T2, x3}, 1588822133503.2803},
      {{T3, x1}, 20278523425.24634},
      {{T3, x2}, 106528081570.3118},
      {{T3, x3}, 1607557369570.622}
    };
    std::map<std::tuple<double, double>, double> ref_KT = {
      {{T1, x1}, 261090124247.34915},
      {{T1, x2}, 583985935874.6938},
      {{T1, x3}, 5143592055967.946},
      {{T2, x1}, 259628514928.99078},
      {{T2, x2}, 583407146677.9891},
      {{T2, x3}, 5142680547136.386},
      {{T3, x1}, 257933159414.78268},
      {{T3, x2}, 586413832596.9236},
      {{T3, x3}, 5152266523051.159}
    };
    std::map<std::tuple<double, double>, double> ref_KT_c = {
      {{T1, x1}, 261090124247.34915},
      {{T1, x2}, 583985935874.6938},
      {{T1, x3}, 5143592055967.946},
      {{T2, x1}, 259596185343.45975},
      {{T2, x2}, 583377935530.4591},
      {{T2, x3}, 5142659545248.286},
      {{T3, x1}, 257571068056.8352},
      {{T3, x2}, 586086667744.5878},
      {{T3, x3}, 5152031301904.44}
    };
    std::map<std::tuple<double, double>, double> ref_KS = {
      {{T1, x1}, 263861771077.48422},
      {{T1, x2}, 585892410050.0793},
      {{T1, x3}, 5144220790689.566},
      {{T2, x1}, 269089945746.04807},
      {{T2, x2}, 591184634672.8638},
      {{T2, x3}, 5149148995109.341},
      {{T3, x1}, 288663292469.2158},
      {{T3, x2}, 612469616901.9039},
      {{T3, x3}, 5179309450284.324}
    };
    std::map<std::tuple<double, double>, double> ref_G2 = {
      {{T1, x1}, 177206739110.75482},
      {{T1, x2}, 297965839231.14703},
      {{T1, x3}, 1582445467995.7842},
      {{T2, x1}, 171415950723.23972},
      {{T2, x2}, 293545321947.4244},
      {{T2, x3}, 1580073907871.653},
      {{T3, x1}, 149736279309.26605},
      {{T3, x2}, 275766289775.94055},
      {{T3, x3}, 1565560034921.9797}
    };
    std::map<std::tuple<double, double>, double> ref_G3 = {
      {{T1, x1}, 177192656138.8814},
      {{T1, x2}, 286530743392.9021},
      {{T1, x3}, 581977273079.2618},
      {{T2, x1}, 171401867751.3663},
      {{T2, x2}, 282110226109.17944},
      {{T2, x3}, 579605712955.1309},
      {{T3, x1}, 149722196337.39264},
      {{T3, x2}, 264331193937.69562},
      {{T3, x3}, 565091840005.4574}
    };
    std::map<std::tuple<double, double>, double> ref_Cv = {
      {{T1, x1}, 93.32288837646946},
      {{T1, x2}, 75.43380616062134},
      {{T1, x3}, 23.317303412816464},
      {{T2, x1}, 119.46436173186248},
      {{T2, x2}, 115.39997421613883},
      {{T2, x3}, 89.95850643101778},
      {{T3, x1}, 124.16407800527138},
      {{T3, x2}, 123.71446776470891},
      {{T3, x3}, 120.34992522915367}
    };
    std::map<std::tuple<double, double>, double> ref_Cv_c = {
      {{T1, x1}, 94.53773268622436},
      {{T1, x2}, 76.32080672751653},
      {{T1, x3}, 23.636166042765097},
      {{T2, x1}, 122.70394655787554},
      {{T2, x2}, 117.76530906119264},
      {{T2, x3}, 90.80880677754747},
      {{T3, x1}, 134.28778058656218},
      {{T3, x2}, 131.10613915550206},
      {{T3, x3}, 123.00711381205896}
    };
    std::map<std::tuple<double, double>, double> ref_Cp = {
      {{T1, x1}, 94.31357344544068},
      {{T1, x2}, 75.6800665490344},
      {{T1, x3}, 23.320153638516846},
      {{T2, x1}, 123.81790430764259},
      {{T2, x2}, 116.93838854511215},
      {{T2, x3}, 90.07165596718679},
      {{T3, x1}, 138.95697491833235},
      {{T3, x2}, 129.21140066141018},
      {{T3, x3}, 120.98161115920557}
    };
    std::map<std::tuple<double, double>, double> ref_S = {
      {{T1, x1}, 70.12285643892892},
      {{T1, x2}, 44.776311261247045},
      {{T1, x3}, 8.634200016536376},
      {{T2, x1}, 177.74371712165933},
      {{T2, x2}, 142.5935489308975},
      {{T2, x3}, 64.2879393541059},
      {{T3, x1}, 317.4605119759535},
      {{T3, x2}, 280.414271100821},
      {{T3, x3}, 189.12990384747835}
    };
    std::map<std::tuple<double, double>, double> ref_S_c = {
      {{T1, x1}, 71.33770074868382},
      {{T1, x2}, 45.66331182814223},
      {{T1, x3}, 8.95306264648501},
      {{T2, x1}, 180.9833019476724},
      {{T2, x2}, 144.9588837759513},
      {{T2, x3}, 65.1382397006356},
      {{T3, x1}, 327.5842145572443},
      {{T3, x2}, 287.80594249161413},
      {{T3, x3}, 191.78709243038364}
    };
    std::map<std::tuple<double, double>, double> ref_F = {
      {{T1, x1}, -1367971.5095417432},
      {{T1, x2}, -1189060.3446439554},
      {{T1, x3}, 3803764.79303171},
      {{T2, x1}, -1433273.5763092036},
      {{T2, x2}, -1238283.8963817155},
      {{T2, x3}, 3786087.8634846513},
      {{T3, x1}, -1875614.8658376308},
      {{T3, x2}, -1618443.3081500707},
      {{T3, x3}, 3555411.020709808}
    };
    std::map<std::tuple<double, double>, double> ref_F_c = {
      {{T1, x1}, -1367971.5095417432},
      {{T1, x2}, -1189060.3446439554},
      {{T1, x3}, 3803764.79303171},
      {{T2, x1}, -1434387.1835931456},
      {{T2, x2}, -1239096.9802347028},
      {{T2, x3}, 3785795.5727405315},
      {{T3, x1}, -1888087.267417781},
      {{T3, x2}, -1627549.8473035279},
      {{T3, x3}, 3552137.3643756686}
    };
    std::map<std::tuple<double, double>, double> ref_E = {
      {{T1, x1}, -1346934.6526100645},
      {{T1, x2}, -1175627.4512655812},
      {{T1, x3}, 3806355.053036671},
      {{T2, x1}, -1291078.602611876},
      {{T2, x2}, -1124209.0572369976},
      {{T2, x3}, 3837518.214967936},
      {{T3, x1}, -1081963.585897747},
      {{T3, x2}, -917407.6303980182},
      {{T3, x3}, 4028235.780328504}
    };
    std::map<std::tuple<double, double>, double> ref_alpha = {
      {{T1, x1}, 2.2861216729857507e-05},
      {{T1, x2}, 8.478039669545874e-06},
      {{T1, x3}, 4.346310245473944e-07},
      {{T2, x1}, 2.9429819049132078e-05},
      {{T2, x2}, 1.298272260723485e-05},
      {{T2, x3}, 1.6771102190256832e-06},
      {{T3, x1}, 3.078863280002657e-05},
      {{T3, x2}, 1.3846757806093788e-05},
      {{T3, x3}, 2.239527608417212e-06}
    };
    std::map<std::tuple<double, double>, double> ref_alpha_c = {
      {{T1, x1}, 2.3145002947829182e-05},
      {{T1, x2}, 8.59267746827458e-06},
      {{T1, x3}, 4.439888174902046e-07},
      {{T2, x1}, 3.0194602491600104e-05},
      {{T2, x2}, 1.3289392083632869e-05},
      {{T2, x3}, 1.7020757075108253e-06},
      {{T3, x1}, 3.322911052986146e-05},
      {{T3, x2}, 1.480637815163179e-05},
      {{T3, x3}, 2.3174837269655006e-06}
    };
    std::map<std::tuple<double, double>, double> ref_gamma_c = {
      {{T1, x1}, 1.5469195895080807},
      {{T1, x2}, 1.2857838750298638},
      {{T1, x3}, 0.9447383595416163},
      {{T2, x1}, 1.5459460283311568},
      {{T2, x2}, 1.287412633212964},
      {{T2, x3}, 0.9425159562230832},
      {{T3, x1}, 1.5424267027049094},
      {{T3, x2}, 1.2943975414783564},
      {{T3, x3}, 0.9491065939171507}
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
        {{P1, T1, x1}, -1409073.0263092036},
        {{P1, T1, x2}, 3795865.8634846513},
        {{P1, T2, x1}, -1851414.3158376308},
        {{P1, T2, x2}, 3565189.020709808},
        {{P2, T1, x1}, -126443.87630920368},
        {{P2, T1, x2}, 4314099.863484651},
        {{P2, T2, x1}, -568785.1658376309},
        {{P2, T2, x2}, 4083423.020709808}
    };
    std::map<std::tuple<double, double, double>, double> ref_H = {
        {{P1, T1, x1}, -1266878.052611876},
        {{P1, T1, x2}, 3847296.214967936},
        {{P1, T2, x1}, -1057763.035897747},
        {{P1, T2, x2}, 4038013.780328504},
        {{P2, T1, x1}, 15751.09738812386},
        {{P2, T1, x2}, 4365530.214967936},
        {{P2, T2, x1}, 224866.114102253},
        {{P2, T2, x2}, 4556247.780328504}
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
