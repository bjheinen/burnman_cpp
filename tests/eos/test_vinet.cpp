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
#include "burnman/eos/vinet.hpp"
#include <cmath>
#include "burnman/utils/eos.hpp"
#include "tolerances.hpp"
using namespace Catch::Matchers;

TEST_CASE("Test validate parameters", "[vinet][eos]") {
  Vinet vinet;
  SECTION("Missing V_0") {
    MineralParams params;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    REQUIRE_THROWS(vinet.validate_parameters(params));
  }
  SECTION("Missing K_0") {
    MineralParams params;
    params.V_0 = 11.24e-6;
    params.Kprime_0 = 3.8;
    REQUIRE_THROWS(vinet.validate_parameters(params));
  }
  SECTION("Missing Kprime_0") {
    MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    REQUIRE_THROWS(vinet.validate_parameters(params));
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
    vinet.validate_parameters(params);
    REQUIRE(*params.E_0 == 0);
    REQUIRE(*params.P_0 == 0);
    REQUIRE(std::isnan(*params.G_0));
    REQUIRE(std::isnan(*params.Gprime_0));
  }
  // TODO: Check warnings
}

TEST_CASE("Check reference volume", "[vinet][eos]") {
  // Set up test params
  MineralParams params;
  params.equation_of_state = EOSType::Vinet;
  params.T_0 = 300.0;
  params.P_0 = 0.0;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.molar_mass = 0.0403;
  params.napfu = 2;
  // Set up test P, T, V
  auto P = GENERATE(0.0, 10.0, 25.e9);
  auto T = GENERATE(300.0, 2000.0);
  double V = *params.V_0;
  Vinet vinet;
  CHECK_THAT(vinet.compute_isothermal_bulk_modulus_reuss(P, T, V, params),
    WithinRel(*params.K_0, tol_rel) || WithinAbs(*params.K_0, tol_abs));
  CHECK_THAT(vinet.compute_isentropic_bulk_modulus_reuss(P, T, V, params),
    WithinRel(*params.K_0, tol_rel) || WithinAbs(*params.K_0, tol_abs));
  CHECK_THAT(vinet.compute_molar_internal_energy(P, T, V, params),
    WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
  CHECK_THAT(vinet.compute_gibbs_free_energy(P, T, V, params),
    WithinRel(P*V, tol_rel) || WithinAbs(P*V, tol_abs));
  CHECK_THAT(vinet.compute_pressure(T, V, params),
    WithinRel(*params.P_0, tol_rel) || WithinAbs(*params.P_0, tol_abs));
  CHECK_THAT(vinet.compute_volume(0.0, T, params),
    WithinRel(V, tol_rel) || WithinAbs(V, tol_abs));
}

TEST_CASE("Check hard-coded returns", "[vinet][eos]") {
  // Set up test P, T, V
  auto P = GENERATE(0.0, 25.e9);
  auto T = GENERATE(300.0, 2000.0);
  auto V = GENERATE(5.0e-6, 10.0e-6);
  // Set up test params
  MineralParams params;
  params.equation_of_state = EOSType::Vinet;
  params.T_0 = 300.0;
  params.P_0 = 0.0;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.molar_mass = 0.0403;
  params.napfu = 2;
  // Unused parameters
  params.G_0 = 131.0e9;
  params.Gprime_0 = 2.1;
  params.grueneisen_0 = 1.5;
  Vinet vinet;
  CHECK(vinet.compute_grueneisen_parameter(P, T, V, params) == 0);
  CHECK(vinet.compute_shear_modulus(P, T, V, params) == 0);
  CHECK(vinet.compute_thermal_expansivity(P, T, V, params) == 0);
  CHECK(vinet.compute_entropy(P, T, V, params) == 0);
  CHECK(vinet.compute_molar_heat_capacity_v(P, T, V, params) == 1.0e99);
  CHECK(vinet.compute_molar_heat_capacity_v(P, T, V, params) == 1.0e99);
}

TEST_CASE("Vinet python reference values", "[vinet][eos]") {
  struct TestData {
    double input;
    double expected_P;
    double expected_K;
    double expected_E;
  };
  // Set up test params
  MineralParams params_a;
  params_a.equation_of_state = EOSType::Vinet;
  params_a.T_0 = 300.0;
  params_a.P_0 = 0.0;
  params_a.E_0 = 0.0;
  params_a.V_0 = 11.24e-6;
  params_a.K_0 = 161.0e9;
  params_a.Kprime_0 = 3.8;
  params_a.molar_mass = 0.0403;
  params_a.napfu = 2;

  MineralParams params_b;
  params_b.equation_of_state = EOSType::Vinet;
  params_b.P_0 = 0.0;
  params_b.E_0 = 0.0;
  params_b.V_0 = 6.75e-6;
  params_b.K_0 = 163.4e9;
  params_b.Kprime_0 = 5.38;
  params_b.molar_mass = 0.055845;
  params_b.napfu = 1;

  SECTION("Test volume dependent functions A") {
    Vinet vinet;
    // P & T unused
    double P = 0.0;
    double T = 300.0;
    auto test_data = GENERATE(
      TestData{0.990, 1649261647.56626, 167226180887.78818, 91.947133283917},
      TestData{0.98, 3379601215.972409, 173674689377.5619, 373.79147656404},
      TestData{0.95, 9098123097.310312, 194449338236.78143, 2454.3925457405157},
      TestData{0.80, 54289565597.9588, 341108449579.4405, 51271.223568070265},
      TestData{0.40, 707279321477.4147, 1861104039241.4866, 1217258.7653399277}
    );
    CAPTURE(test_data.input);
    CHECK_THAT(vinet.compute_pressure(
        T, test_data.input*(*params_a.V_0), params_a
      ),
      WithinRel(test_data.expected_P, tol_rel) ||
      WithinAbs(test_data.expected_P, tol_abs));
    CHECK_THAT(vinet.compute_isothermal_bulk_modulus_reuss(
        P, T, test_data.input*(*params_a.V_0), params_a
      ),
      WithinRel(test_data.expected_K, tol_rel) ||
      WithinAbs(test_data.expected_K, tol_abs));
    CHECK_THAT(vinet.compute_isentropic_bulk_modulus_reuss(
        P, T, test_data.input*(*params_a.V_0), params_a
      ),
      WithinRel(test_data.expected_K, tol_rel) ||
      WithinAbs(test_data.expected_K, tol_abs));
    CHECK_THAT(vinet.compute_molar_internal_energy(
        P, T, test_data.input*(*params_a.V_0), params_a
      ),
      WithinRel(test_data.expected_E, tol_rel) ||
      WithinAbs(test_data.expected_E, tol_abs));
  }
  SECTION("Test volume dependent functions B") {
    Vinet vinet;
    // P & T unused
    double P = 20.e9;
    double T = 2000.0;
    auto test_data = GENERATE(
      TestData{0.99, 1687167351.8287652, 172398018035.712, 56.33784395172695},
      TestData{0.98, 3484975646.432192, 181824440103.00647, 230.2567007936093},
      TestData{0.95, 9612287325.901537, 212903271610.8321, 1536.838065740325},
      TestData{0.80, 65301548960.54173, 458188447601.4275, 35120.785814691655},
      TestData{0.40, 1339419452737.9487, 4304135307666.401, 1175069.4080477764}
    );
    CAPTURE(test_data.input);
    CHECK_THAT(vinet.compute_pressure(
        T, test_data.input*(*params_b.V_0), params_b
      ),
      WithinRel(test_data.expected_P, tol_rel) ||
      WithinAbs(test_data.expected_P, tol_abs));
    CHECK_THAT(vinet.compute_isothermal_bulk_modulus_reuss(
        P, T, test_data.input*(*params_b.V_0), params_b
      ),
      WithinRel(test_data.expected_K, tol_rel) ||
      WithinAbs(test_data.expected_K, tol_abs));
    CHECK_THAT(vinet.compute_isentropic_bulk_modulus_reuss(
        P, T, test_data.input*(*params_b.V_0), params_b
      ),
      WithinRel(test_data.expected_K, tol_rel) ||
      WithinAbs(test_data.expected_K, tol_abs));
    CHECK_THAT(vinet.compute_molar_internal_energy(
        P, T, test_data.input*(*params_b.V_0), params_b
      ),
      WithinRel(test_data.expected_E, tol_rel) ||
      WithinAbs(test_data.expected_E, tol_abs));
  }
  SECTION("Test Gibbs") {
    Vinet vinet;
    double T = 2000.0; // Unused
    double V_aa = *params_a.V_0 * 0.99;
    double P_aa = 1.6e9;
    double V_ab = *params_b.V_0 * 0.99;
    double P_ab = 1.6e9;
    double V_ba = *params_a.V_0 * 0.8;
    double P_ba = 54.0e9;
    double V_bb = *params_b.V_0 * 0.8;
    double P_bb = 65.0e9;
    double expected_G_aa = 17896.10713328402;
    double expected_G_ab = 10748.337843951727;
    double expected_G_ba = 536839.2235680703;
    double expected_G_bb = 386120.7858146917;
    CHECK_THAT(vinet.compute_gibbs_free_energy(P_aa, T, V_aa, params_a),
      WithinRel(expected_G_aa, tol_rel) || WithinAbs(expected_G_aa, tol_abs));
    CHECK_THAT(vinet.compute_gibbs_free_energy(P_ab, T, V_ab, params_b),
      WithinRel(expected_G_ab, tol_rel) || WithinAbs(expected_G_ab, tol_abs));
    CHECK_THAT(vinet.compute_gibbs_free_energy(P_ba, T, V_ba, params_a),
      WithinRel(expected_G_ba, tol_rel) || WithinAbs(expected_G_ba, tol_abs));
    CHECK_THAT(vinet.compute_gibbs_free_energy(P_bb, T, V_bb, params_b),
      WithinRel(expected_G_bb, tol_rel) || WithinAbs(expected_G_bb, tol_abs));
  }
  SECTION("Test volume") {
    Vinet vinet;
    double T = 2000.0;
    double V_a = 0.9 * (*params_a.V_0);
    double V_b = 0.5 * (*params_b.V_0);
    double P_a = vinet.compute_pressure(T, V_a, params_a);
    double P_b = vinet.compute_pressure(T, V_b, params_b);
    CHECK_THAT(vinet.compute_volume(P_a, T, params_a),
      WithinRel(V_a, tol_rel) || WithinAbs(V_a, tol_abs));
    CHECK_THAT(vinet.compute_volume(P_b, T, params_b),
      WithinRel(V_b, tol_rel) || WithinAbs(V_b, tol_abs));
  }
}
