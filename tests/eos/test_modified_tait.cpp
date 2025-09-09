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
#include <cmath>
#include <map>
#include <tuple>
#include "burnman/utils/eos.hpp"
#include "burnman/eos/modified_tait.hpp"
#include "tolerances.hpp"
using namespace Catch::Matchers;

TEST_CASE("Test validate parameters", "[modified_tait][eos]") {
  MineralParams params;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.Kdprime_0 = -4e-11;
  params.G_0 = 125.0e9;
  params.Gprime_0 = 1.8;
  params.E_0 = 0.0;
  params.P_0 = 1.0e5;
  MT mt;
  SECTION("Missing V_0") {
    params.V_0.reset();
    REQUIRE_THROWS(mt.validate_parameters(params));
  }
  SECTION("Missing K_0") {
    params.K_0.reset();
    REQUIRE_THROWS(mt.validate_parameters(params));
  }
  SECTION("Missing Kprime_0") {
    params.Kprime_0.reset();
    REQUIRE_THROWS(mt.validate_parameters(params));
  }
  SECTION("Missing Kdprime_0") {
    params.Kdprime_0.reset();
    REQUIRE_THROWS(mt.validate_parameters(params));
  }
  SECTION("Optional Parameters") {
    params.E_0.reset();
    params.P_0.reset();
    params.G_0.reset();
    params.Gprime_0.reset();
    REQUIRE_FALSE(params.E_0.has_value());
    REQUIRE_FALSE(params.P_0.has_value());
    REQUIRE_FALSE(params.G_0.has_value());
    REQUIRE_FALSE(params.Gprime_0.has_value());
    REQUIRE_NOTHROW(mt.validate_parameters(params));
    REQUIRE(*params.E_0 == 0);
    REQUIRE(*params.P_0 == 1.0e5);
    REQUIRE(std::isnan(*params.G_0));
    REQUIRE(std::isnan(*params.Gprime_0));
  }
  // TODO: Check warnings
}

TEST_CASE("Tait constant", "[modified_tait][eos]") {
  MineralParams params;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 4.0;
  params.Kdprime_0 = -4.0e-11;
  TaitConstants result;
  SECTION("Kdprime_0 zero behaviour") {
    params.Kdprime_0 = 0;
    result = MT::compute_tait_constants(params);
    CHECK(result.a == 1.0);
  }
  SECTION("Kprime_0 zero behaviour") {
    params.Kprime_0 = 0;
    result = MT::compute_tait_constants(params);
    CHECK_THAT(result.b,
      WithinRel(-*params.Kdprime_0, tol_rel) ||
      WithinAbs(-*params.Kdprime_0, tol_abs));
  }
  SECTION("Reference check") {
    double ref_a = -625.0/180.0;
    double ref_b = 3.284472049689441e-11;
    double ref_c = -0.05446293494704991;
    result = MT::compute_tait_constants(params);
    CHECK_THAT(result.a,
      WithinRel(ref_a, tol_rel) ||
      WithinAbs(ref_a, tol_abs));
    CHECK_THAT(result.b,
      WithinRel(ref_b, tol_rel) ||
      WithinAbs(ref_b, tol_abs));
    CHECK_THAT(result.c,
      WithinRel(ref_c, tol_rel) ||
      WithinAbs(ref_c, tol_abs));
  }
}

TEST_CASE("Check reference conditions", "[modified_tait][eos]") {
  // Set up test params
  MineralParams params;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.Kdprime_0 = -4e-11;
  params.G_0 = 125.0e9;
  params.Gprime_0 = 1.8;
  params.E_0 = 0.0;
  params.P_0 = 1.0e5;
  MT mt;
  SECTION("Reference V") {
    auto T = GENERATE(300.0, 2000.0);
    double V = *params.V_0;
    CHECK_THAT(mt.compute_pressure(T, V, params),
      WithinRel(*params.P_0, tol_rel) || WithinAbs(*params.P_0, tol_abs));
  }
  SECTION("Reference P") {
    auto P = *params.P_0;
    auto T = GENERATE(300.0, 2000.0);
    auto V = GENERATE(11.24e-6, 11.24e-6*0.8, 11.24e-6*0.4);
    CHECK_THAT(mt.compute_isothermal_bulk_modulus_reuss(P, T, V, params),
      WithinRel(*params.K_0, tol_rel) || WithinAbs(*params.K_0, tol_abs));
    CHECK_THAT(mt.compute_gibbs_free_energy(P, T, V, params),
      WithinRel((*params.P_0)*(*params.V_0), tol_rel)
        || WithinAbs((*params.P_0)*(*params.V_0), tol_abs));
    CHECK_THAT(mt.compute_volume(*params.P_0, T, params),
      WithinRel(*params.V_0, tol_rel) || WithinAbs(*params.V_0, tol_abs));
    // K_S returns 1.0e99 in MT
    //CHECK_THAT(mt.compute_isentropic_bulk_modulus_reuss(P, T, V, params),
    //  WithinRel(*params.K_0, tol_rel) || WithinAbs(*params.K_0, tol_abs));
  }
  SECTION("Reference P-V") {
    auto T = GENERATE(300.0, 2000.0);
    double P = *params.P_0;
    double V = *params.V_0;
    CHECK_THAT(mt.compute_molar_internal_energy(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
  }
}

TEST_CASE("Check hard-coded returns", "[modified_tait][eos]") {
  // Set up test P, T, V
  auto P = GENERATE(0.0, 25.e9);
  auto T = GENERATE(300.0, 2000.0);
  auto V = GENERATE(5.0e-6, 10.0e-6);
  // Set up test params
  MineralParams params;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.Kdprime_0 = -4e-11;
  params.G_0 = 125.0e9;
  params.Gprime_0 = 1.8;
  params.E_0 = 0.0;
  params.P_0 = 1.0e5;
  MT mt;
  CHECK(mt.compute_thermal_expansivity(P, T, V, params) == 0);
  CHECK(mt.compute_grueneisen_parameter(P, T, V, params) == 0);
  CHECK(mt.compute_molar_heat_capacity_v(P, T, V, params) == 1.0e99);
  CHECK(mt.compute_molar_heat_capacity_v(P, T, V, params) == 1.0e99);
  CHECK(mt.compute_entropy(P, T, V, params) == 0);
  CHECK(mt.compute_shear_modulus(P, T, V, params) == 0);
  CHECK(mt.compute_isentropic_bulk_modulus_reuss(P, T, V, params) == 1.0e99);
}

TEST_CASE("MT python reference values", "[modified_tait][eos]") {
  // Set up test params
  MineralParams params;
  params.V_0 = 2.445e-05;
  params.K_0 = 251.0e9;
  params.Kprime_0 = 4.14;
  params.Kdprime_0 = -1.6e-11;
  MT mt;
  mt.validate_parameters(params);
  SECTION("Test volume dependent functions") {
    struct TestData {
      double input;
      double expected_P;
    };
    double T = 300.0;
    auto test_data = GENERATE(
      TestData{0.99, 2575775494.857839},
      TestData{0.98, 5287609064.26949},
      TestData{0.95, 14317458932.430204},
      TestData{0.80, 88532160429.76569},
      TestData{0.40, 1157563719136.7415}
    );
    CHECK_THAT(mt.compute_pressure(
        T, test_data.input*(*params.V_0), params
      ),
      WithinRel(test_data.expected_P, tol_rel) ||
      WithinAbs(test_data.expected_P, tol_abs));
  }
  SECTION("Test pressure dependent functions") {
    struct TestData {
      double input;
      double expected_K;
      double expected_G;
    };
    // P & T unused
    double V = *params.V_0;
    double T = 300.0;
    auto test_data = GENERATE(
      TestData{0.0, 250999585999.92, -4.868470204755226e-07},
      TestData{1.0e9, 255131637073.39673, 24401.633633392434},
      TestData{5.5e9, 273535494608.45435, 133054.17058804224},
      TestData{25.0e9, 350134349013.2772, 584996.0112561124},
      TestData{125.0e9, 691061592418.0354, 2603347.111664467}
    );
    CAPTURE(test_data.input);
    CHECK_THAT(mt.compute_isothermal_bulk_modulus_reuss(
        test_data.input, T, V, params
      ),
      WithinRel(test_data.expected_K, tol_rel) ||
      WithinAbs(test_data.expected_K, tol_abs));
    CHECK_THAT(mt.compute_gibbs_free_energy(
        test_data.input, T, V, params
      ),
      WithinRel(test_data.expected_G, tol_rel) ||
      WithinAbs(test_data.expected_G, tol_abs));
  }
  SECTION("Test P-V depedent") {
    double T = 300.0;
    double P1 = 1.0e9, P2 = 25.0e9, P3 = 125.0e9;
    double x1 = 0.99, x2 = 0.80, x3 = 0.40;
    // Reference data
    std::map<std::tuple<double, double>, double> ref_E = {
      {{P1, x1}, 196.1336333924337},
      {{P1, x2}, 4841.63363339243},
      {{P1, x3}, 14621.633633392432},
      {{P2, x1}, -20141.488743887632},
      {{P2, x2}, 95996.01125611231},
      {{P2, x3}, 340496.01125611237},
      {{P3, x1}, -422340.388335533},
      {{P3, x2}, 158347.11166446656},
      {{P3, x3}, 1380847.1116644668}
    };
    auto P = GENERATE(1.0e9, 25.0e9, 125.0e9);
    auto x = GENERATE(0.99, 0.8, 0.4);
    CAPTURE(P);
    CAPTURE(x);
    double V = *params.V_0 * x;
    auto key = std::make_tuple(P, x);
    CHECK_THAT(mt.compute_molar_internal_energy(P, T, V, params),
      WithinRel(ref_E[key], tol_rel) || WithinAbs(ref_E[key], tol_abs));
  }
  SECTION("Test volume") {
    double T = 2000.0;
    double V_a = 0.9 * (*params.V_0);
    double V_b = 0.4 * (*params.V_0);
    double P_a = mt.compute_pressure(T, V_a, params);
    double P_b = mt.compute_pressure(T, V_b, params);
    CHECK_THAT(mt.compute_volume(P_a, T, params),
      WithinRel(V_a, tol_rel) || WithinAbs(V_a, tol_abs));
    CHECK_THAT(mt.compute_volume(P_b, T, params),
      WithinRel(V_b, tol_rel) || WithinAbs(V_b, tol_abs));
  }
}
