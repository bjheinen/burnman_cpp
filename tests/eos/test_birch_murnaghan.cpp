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
#include "burnman/eos/birch_murnaghan.hpp"
#include <cmath>
#include "burnman/utils/types/mineral_params.hpp"
#include "tolerances.hpp"
using namespace Catch::Matchers;

TEST_CASE("Test validate parameters", "[birch-murnaghan][eos]") {
  BM3 bm3;
  SECTION("Missing V_0") {
    MineralParams params;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    REQUIRE_THROWS(bm3.validate_parameters(params));
  }
  SECTION("Missing K_0") {
    MineralParams params;
    params.V_0 = 11.24e-6;
    params.Kprime_0 = 3.8;
    REQUIRE_THROWS(bm3.validate_parameters(params));
  }
  SECTION("Missing Kprime_0") {
    MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    REQUIRE_THROWS(bm3.validate_parameters(params));
  }
  SECTION("Optional Parameters") {
    MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    REQUIRE_FALSE(params.E_0.has_value());
    REQUIRE_FALSE(params.P_0.has_value());
    REQUIRE_FALSE(params.G_0.has_value());
    REQUIRE_FALSE(params.Gprime_0.has_value());
    bm3.validate_parameters(params);
    REQUIRE(*params.E_0 == 0);
    REQUIRE(*params.P_0 == 0);
    REQUIRE(std::isnan(*params.G_0));
    REQUIRE(std::isnan(*params.Gprime_0));
  }
  // TODO: Check warnings
}

TEST_CASE("Check reference volume", "[birch-murnaghan][eos]") {
  // Set up test params
  MineralParams params_a;
  params_a.equation_of_state = EOSType::BM3;
  params_a.T_0 = 300.0;
  params_a.P_0 = 0.0;
  params_a.V_0 = 11.24e-6;
  params_a.K_0 = 161.0e9;
  params_a.Kprime_0 = 3.8;
  params_a.molar_mass = 0.0403;
  params_a.napfu = 2;
  MineralParams params_b = params_a;
  params_a.G_0 = 131.0e9;
  params_b.G_0 = 129.8e9;
  params_a.Gprime_0 = 2.1;
  // Set up test P, T, V
  auto P = GENERATE(0.0, 10.0, 25.e9);
  auto T = GENERATE(300.0, 2000.0);
  double V = *params_a.V_0;
  BM3 bm3;
  BM2 bm2;
  CHECK_THAT(bm3.compute_isothermal_bulk_modulus_reuss(P, T, V, params_a),
    WithinRel(*params_a.K_0, tol_rel) || WithinAbs(*params_a.K_0, tol_abs));
  CHECK_THAT(bm3.compute_isentropic_bulk_modulus_reuss(P, T, V, params_a),
    WithinRel(*params_a.K_0, tol_rel) || WithinAbs(*params_a.K_0, tol_abs));
  CHECK_THAT(bm3.compute_molar_internal_energy(P, T, V, params_a),
    WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
  CHECK_THAT(bm3.compute_gibbs_free_energy(P, T, V, params_a),
    WithinRel(P*V, tol_rel) || WithinAbs(P*V, tol_abs));
  CHECK_THAT(bm3.compute_pressure(T, V, params_a),
    WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
  CHECK_THAT(bm3.compute_volume(0.0, T, params_a),
    WithinRel(V, tol_rel) || WithinAbs(V, tol_abs));
  // Check BM3 and BM2 shear modulus
  CHECK_THAT(bm3.compute_shear_modulus(P, T, V, params_a),
    WithinRel(*params_a.G_0, tol_rel) || WithinAbs(*params_a.G_0, tol_abs));
  CHECK_THAT(bm2.compute_shear_modulus(P, T, V, params_b),
    WithinRel(*params_b.G_0, tol_rel) || WithinAbs(*params_b.G_0, tol_abs));
}

TEST_CASE("Check hard-coded returns", "[birch-murnaghan][eos]") {
  // Set up test P, T, V
  auto P = GENERATE(0.0, 25.e9);
  auto T = GENERATE(300.0, 2000.0);
  auto V = GENERATE(5.0e-6, 10.0e-6);
  // Set up test params
  MineralParams params;
  params.equation_of_state = EOSType::BM2;
  params.T_0 = 300.0;
  params.P_0 = 0.0;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.molar_mass = 0.0403;
  params.napfu = 2;
  // Unused parameters
  params.G_0 = std::nan("");
  params.Gprime_0 = std::nan("");
  params.grueneisen_0 = 1.5;
  BM2 bm2;
  CHECK(bm2.compute_grueneisen_parameter(P, T, V, params) == 0);
  CHECK(bm2.compute_thermal_expansivity(P, T, V, params) == 0);
  CHECK(bm2.compute_entropy(P, T, V, params) == 0);
  CHECK(bm2.compute_molar_heat_capacity_v(P, T, V, params) == 1.0e99);
  CHECK(bm2.compute_molar_heat_capacity_v(P, T, V, params) == 1.0e99);
}

TEST_CASE("Shear modulus expansion", "[birch-murnaghan][eos]") {
  // Set up test P, T, V
  double P = 25.e9;
  double T = 2000.0;
  double V = 10.0e-6;
  // Set up test params
  MineralParams params;
  params.T_0 = 300.0;
  params.P_0 = 0.0;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.molar_mass = 0.0403;
  params.napfu = 2;
  BM2 bm2;
  BM3 bm3;
  // Validating should set G0/G' to nan
  bm3.validate_parameters(params);
  REQUIRE_NOTHROW(bm2.compute_shear_modulus(P, T, V, params));
  REQUIRE_NOTHROW(bm3.compute_shear_modulus(P, T, V, params));
  CHECK(std::isnan(bm2.compute_shear_modulus(P, T, V, params)));
  CHECK(std::isnan(bm3.compute_shear_modulus(P, T, V, params)));
  // Set G0/G'0
  params.G_0 = 131.0e9;
  params.Gprime_0 = 2.1;
  CHECK_FALSE(std::isnan(bm3.compute_shear_modulus(P, T, V, params)));
  // Ensure expansion order different
  CHECK_FALSE(
    bm2.compute_shear_modulus(P, T, V, params) ==
    bm3.compute_shear_modulus(P, T, V, params)
  );
}

TEST_CASE("Birch-Murnaghan python reference values", "[birch-murnaghan][eos]") {
  struct TestData {
    double input;
    double expected_P;
    double expected_K;
    double expected_E;
    double expected_G2;
    double expected_G3;
  };
  // Set up test params
  MineralParams params_a;
  params_a.T_0 = 300.0;
  params_a.P_0 = 0.0;
  params_a.E_0 = 0.0;
  params_a.V_0 = 11.24e-6;
  params_a.K_0 = 161.0e9;
  params_a.Kprime_0 = 3.8;
  params_a.molar_mass = 0.0403;
  params_a.napfu = 2;
  params_a.G_0 = 131.0e9;
  params_a.Gprime_0 = 2.1;

  MineralParams params_b;
  params_b.P_0 = 0.0;
  params_b.E_0 = 0.0;
  params_b.V_0 = 6.75e-6;
  params_b.K_0 = 163.4e9;
  params_b.Kprime_0 = 5.38;
  params_b.molar_mass = 0.055845;
  params_b.napfu = 1;
  params_b.G_0 = 124.3e9;
  params_b.Gprime_0 = 1.7;

  SECTION("Test volume dependent functions A") {
    BM3 bm3;
    BM2 bm2;
    // P & T unused
    double P = 0.0;
    double T = 300.0;
    auto test_data = GENERATE(
      TestData{0.99, 1649296384.9124484, 167236609646.01492, 91.94810310156963, 134440941855.59421, 134430413937.49939},
      TestData{0.98, 3379888334.227698, 173717817989.04596, 373.80740393158857, 138005422822.76697, 137961862270.56598},
      TestData{0.95, 9103079175.634022, 194747670891.14215, 2455.066251924619, 149498901056.86096, 149196949733.73715},
      TestData{0.80, 54834505788.6455, 349338640343.3208, 51536.61471078358, 231811619814.15884, 223263157637.72803},
      TestData{0.40, 818147973581.2825, 2384291228322.9097, 1321858.6932650988, 1299855849555.416, 551942177125.8463}
    );
    CAPTURE(test_data.input);
    CHECK_THAT(bm3.compute_pressure(
        T, test_data.input*(*params_a.V_0), params_a
      ),
      WithinRel(test_data.expected_P, tol_rel) ||
      WithinAbs(test_data.expected_P, tol_abs));
    CHECK_THAT(bm3.compute_isothermal_bulk_modulus_reuss(
        P, T, test_data.input*(*params_a.V_0), params_a
      ),
      WithinRel(test_data.expected_K, tol_rel) ||
      WithinAbs(test_data.expected_K, tol_abs));
    CHECK_THAT(bm3.compute_isentropic_bulk_modulus_reuss(
        P, T, test_data.input*(*params_a.V_0), params_a
      ),
      WithinRel(test_data.expected_K, tol_rel) ||
      WithinAbs(test_data.expected_K, tol_abs));
    CHECK_THAT(bm3.compute_molar_internal_energy(
        P, T, test_data.input*(*params_a.V_0), params_a
      ),
      WithinRel(test_data.expected_E, tol_rel) ||
      WithinAbs(test_data.expected_E, tol_abs));
    CHECK_THAT(bm2.compute_shear_modulus(
        P, T, test_data.input*(*params_a.V_0), params_a
      ),
      WithinRel(test_data.expected_G2, tol_rel) ||
      WithinAbs(test_data.expected_G2, tol_abs));
    CHECK_THAT(bm3.compute_shear_modulus(
        P, T, test_data.input*(*params_a.V_0), params_a
      ),
      WithinRel(test_data.expected_G3, tol_rel) ||
      WithinAbs(test_data.expected_G3, tol_abs));
  }
  SECTION("Test volume dependent functions B") {
    BM3 bm3;
    BM2 bm2;
    // P & T unused
    double P = 20.e9;
    double T = 2000.0;
    auto test_data = GENERATE(
      TestData{0.99, 1687230485.0584817, 172417026177.66202, 56.33890061649582, 127123729083.17744, 127123278622.02191},
      TestData{0.98, 3485502055.263302, 181903966690.09323, 230.27417654703999, 130042001212.17201, 130040137373.89894},
      TestData{0.95, 9621620420.08344, 213473107960.34708, 1537.5933163992263, 139407180309.35876, 139394260628.49622},
      TestData{0.80, 66490483921.849304, 477159619815.7431, 35455.90643881672, 204939596463.2018, 204573830875.41495},
      TestData{0.40, 1778621452528.4504, 6766264877390.478, 1390843.0929669242, 983112854774.1771, 951111650538.3658}
    );
    CAPTURE(test_data.input);
    CHECK_THAT(bm3.compute_pressure(
        T, test_data.input*(*params_b.V_0), params_b
      ),
      WithinRel(test_data.expected_P, tol_rel) ||
      WithinAbs(test_data.expected_P, tol_abs));
    CHECK_THAT(bm3.compute_isothermal_bulk_modulus_reuss(
        P, T, test_data.input*(*params_b.V_0), params_b
      ),
      WithinRel(test_data.expected_K, tol_rel) ||
      WithinAbs(test_data.expected_K, tol_abs));
    CHECK_THAT(bm3.compute_isentropic_bulk_modulus_reuss(
        P, T, test_data.input*(*params_b.V_0), params_b
      ),
      WithinRel(test_data.expected_K, tol_rel) ||
      WithinAbs(test_data.expected_K, tol_abs));
    CHECK_THAT(bm3.compute_molar_internal_energy(
        P, T, test_data.input*(*params_b.V_0), params_b
      ),
      WithinRel(test_data.expected_E, tol_rel) ||
      WithinAbs(test_data.expected_E, tol_abs));
    CHECK_THAT(bm2.compute_shear_modulus(
        P, T, test_data.input*(*params_b.V_0), params_b
      ),
      WithinRel(test_data.expected_G2, tol_rel) ||
      WithinAbs(test_data.expected_G2, tol_abs));
    CHECK_THAT(bm3.compute_shear_modulus(
        P, T, test_data.input*(*params_b.V_0), params_b
      ),
      WithinRel(test_data.expected_G3, tol_rel) ||
      WithinAbs(test_data.expected_G3, tol_abs));
  }
  SECTION("Test Gibbs") {
    BM3 bm3;
    double T = 2000.0; // Unused
    double V_aa = *params_a.V_0 * 0.99;
    double P_aa = 1.6e9;
    double V_ab = *params_b.V_0 * 0.99;
    double P_ab = 1.6e9;
    double V_ba = *params_a.V_0 * 0.8;
    double P_ba = 54.0e9;
    double V_bb = *params_b.V_0 * 0.8;
    double P_bb = 65.0e9;
    double expected_G_aa = 17896.10810310157;
    double expected_G_ab = 10748.338900616496;
    double expected_G_ba = 537104.6147107836;
    double expected_G_bb = 386455.9064388167;
    CHECK_THAT(bm3.compute_gibbs_free_energy(P_aa, T, V_aa, params_a),
      WithinRel(expected_G_aa, tol_rel) || WithinAbs(expected_G_aa, tol_abs));
    CHECK_THAT(bm3.compute_gibbs_free_energy(P_ab, T, V_ab, params_b),
      WithinRel(expected_G_ab, tol_rel) || WithinAbs(expected_G_ab, tol_abs));
    CHECK_THAT(bm3.compute_gibbs_free_energy(P_ba, T, V_ba, params_a),
      WithinRel(expected_G_ba, tol_rel) || WithinAbs(expected_G_ba, tol_abs));
    CHECK_THAT(bm3.compute_gibbs_free_energy(P_bb, T, V_bb, params_b),
      WithinRel(expected_G_bb, tol_rel) || WithinAbs(expected_G_bb, tol_abs));
  }
  SECTION("Test volume") {
    BM2 bm2;
    double T = 2000.0;
    double V_a = 0.9 * (*params_a.V_0);
    double V_b = 0.5 * (*params_b.V_0);
    double P_a = bm2.compute_pressure(T, V_a, params_a);
    double P_b = bm2.compute_pressure(T, V_b, params_b);
    CHECK_THAT(bm2.compute_volume(P_a, T, params_a),
      WithinRel(V_a, tol_rel) || WithinAbs(V_a, tol_abs));
    CHECK_THAT(bm2.compute_volume(P_b, T, params_b),
      WithinRel(V_b, tol_rel) || WithinAbs(V_b, tol_abs));
  }
}
