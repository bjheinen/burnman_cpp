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
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "burnman/eos/debye.hpp"
#include <cmath>
#include "tolerances.hpp"

using namespace debye;
using namespace Catch::Matchers;

TEST_CASE("Debye function integrand", "[debye][eos]") {
  REQUIRE(std::isnan(debye_fn_integrand(0.0, nullptr)));
  double xi = std::log(2);
  double ref_log2 = xi * xi * xi;
  CHECK_THAT(debye_fn_integrand(xi, nullptr),
    WithinRel(ref_log2, tol_rel) || WithinAbs(ref_log2, tol_abs));
  xi = 2;
  double ref_2 = 8 / (std::exp(2) - 1);
  CHECK_THAT(debye_fn_integrand(xi, nullptr),
    WithinRel(ref_2, tol_rel) || WithinAbs(ref_2, tol_abs));
}

TEST_CASE("Debye function", "[debye][eos]") {
  // Test values from GSL test suite
  double ref_a = 0.962999940487211048;
  double ref_b = 0.674415564077814667;
  double ref_c = 0.0192957656903454886;
  SECTION("GSL") {
    REQUIRE_THAT(debye_fn_cheb(0.1),
      WithinRel(ref_a, tol_rel) || WithinAbs(ref_a, tol_abs));
    REQUIRE_THAT(debye_fn_cheb(1.0),
      WithinRel(ref_b, tol_rel) || WithinAbs(ref_b, tol_abs));
    REQUIRE_THAT(debye_fn_cheb(10.0),
      WithinRel(ref_c, tol_rel) || WithinAbs(ref_c, tol_abs));
  }
  SECTION("QUAD") {
    CHECK_THAT(debye_fn_quad(0.1),
      WithinRel(ref_a, tol_rel) || WithinAbs(ref_a, tol_abs));
    CHECK_THAT(debye_fn_quad(1.0),
      WithinRel(ref_b, tol_rel) || WithinAbs(ref_b, tol_abs));
    CHECK_THAT(debye_fn_quad(10.0),
      WithinRel(ref_c, tol_rel) || WithinAbs(ref_c, tol_abs));
  }
  SECTION("Debye fn equivalence") {
    double quad_val = debye_fn_quad(300.0);
    CHECK_THAT(debye_fn_cheb(300.0),
      WithinRel(quad_val, tol_rel) || WithinAbs(quad_val, tol_abs));
  }
}

TEST_CASE("ExplicitDouble overloads", "[eos][debye]") {
  double theta_0 = 773.0;
  double T = 1000.0;
  int int_napfu = 2;
  ExplicitDouble dbl_napfu = ExplicitDouble(2.0);
  // Using int
  double int_E = compute_thermal_energy(T, theta_0, int_napfu);
  double int_C = compute_molar_heat_capacity_v(T, theta_0, int_napfu);
  double int_F = compute_helmholtz_free_energy(T, theta_0, int_napfu);
  double int_S = compute_entropy(T, theta_0, int_napfu);
  double int_dCdT = compute_dmolar_heat_capacity_v_dT(T, theta_0, int_napfu);
  // Using ExplicitDouble
  double dbl_E = compute_thermal_energy(T, theta_0, dbl_napfu);
  double dbl_C = compute_molar_heat_capacity_v(T, theta_0, dbl_napfu);
  double dbl_F = compute_helmholtz_free_energy(T, theta_0, dbl_napfu);
  double dbl_S = compute_entropy(T, theta_0, dbl_napfu);
  double dbl_dCdT = compute_dmolar_heat_capacity_v_dT(T, theta_0, dbl_napfu);
  CHECK_THAT(int_E,
    WithinRel(dbl_E, tol_rel) || WithinAbs(dbl_E, tol_abs));
  CHECK_THAT(int_C,
    WithinRel(dbl_C, tol_rel) || WithinAbs(dbl_C, tol_abs));
  CHECK_THAT(int_F,
    WithinRel(dbl_F, tol_rel) || WithinAbs(dbl_F, tol_abs));
  CHECK_THAT(int_S,
    WithinRel(dbl_S, tol_rel) || WithinAbs(dbl_S, tol_abs));
  CHECK_THAT(int_dCdT,
    WithinRel(dbl_dCdT, tol_rel) || WithinAbs(dbl_dCdT, tol_abs));
}

TEST_CASE("Check zero returns in debye model functions", "[debye][eos]") {
  double debye_0 = 773.0;
  int napfu = 2;
  // Explicit return of 0
  CHECK(compute_thermal_energy(0, debye_0, napfu) == 0);
  CHECK(compute_molar_heat_capacity_v(0, debye_0, napfu) == 0);
  CHECK(compute_helmholtz_free_energy(0, debye_0, napfu) == 0);
  CHECK(compute_entropy(0, debye_0, napfu) == 0);
  CHECK(compute_dmolar_heat_capacity_v_dT(0, debye_0, napfu) == 0);

  double small_val = 1.e-16;
  CHECK(compute_thermal_energy(small_val, debye_0, napfu) == 0);
  CHECK(compute_molar_heat_capacity_v(small_val, debye_0, napfu) == 0);
  CHECK(compute_helmholtz_free_energy(small_val, debye_0, napfu) == 0);
  CHECK(compute_entropy(small_val, debye_0, napfu) == 0);
  CHECK(compute_dmolar_heat_capacity_v_dT(small_val, debye_0, napfu) == 0);
}

// Reference values from Py burnman v2.1.1a0
TEST_CASE("Debye napfu constant", "[debye][eos]") {
  double T = 800.0;
  double debye = 500.0;
  double check_half_n;
  check_half_n = compute_thermal_energy(T, debye, 1);
  CHECK_THAT(compute_thermal_energy(T, debye, 2) / 2.0,
    WithinRel(check_half_n, tol_rel) || WithinAbs(check_half_n, tol_abs));
  check_half_n = compute_molar_heat_capacity_v(T, debye, 1);
  CHECK_THAT(compute_molar_heat_capacity_v(T, debye, 2) / 2.0,
    WithinRel(check_half_n, tol_rel) || WithinAbs(check_half_n, tol_abs));
  check_half_n = compute_helmholtz_free_energy(T, debye, 1);
  CHECK_THAT(compute_helmholtz_free_energy(T, debye, 2) / 2.0,
    WithinRel(check_half_n, tol_rel) || WithinAbs(check_half_n, tol_abs));
  check_half_n = compute_entropy(T, debye, 1);
  CHECK_THAT(compute_entropy(T, debye, 2) / 2.0,
    WithinRel(check_half_n, tol_rel) || WithinAbs(check_half_n, tol_abs));
  check_half_n = compute_dmolar_heat_capacity_v_dT(T, debye, 1);
  CHECK_THAT(compute_dmolar_heat_capacity_v_dT(T, debye, 2) / 2.0,
    WithinRel(check_half_n, tol_rel) || WithinAbs(check_half_n, tol_abs));
}

TEST_CASE("Debye functions python reference values", "[debye][eos]") {
  double T_a = 2000.0;
  double debye_a = 543.0;
  int n_a = 2;
  double T_b = 300.0;
  double debye_b = 773.0;
  int n_b = 1;
  SECTION("compute_thermal_energy") {
    double ref_a = 89982.76111222361;
    double ref_b = 2561.5489877441814;
    CHECK_THAT(compute_thermal_energy(T_a, debye_a, n_a),
      WithinRel(ref_a, tol_rel) || WithinAbs(ref_a, tol_abs));
    CHECK_THAT(compute_thermal_energy(T_b, debye_b, n_b),
      WithinRel(ref_b, tol_rel) || WithinAbs(ref_b, tol_abs));
  }
  SECTION("compute_molar_heat_capacity_v") {
    double ref_a = 49.703395320864054;
    double ref_b = 18.288859954203573;
    CHECK_THAT(compute_molar_heat_capacity_v(T_a, debye_a, n_a),
      WithinRel(ref_a, tol_rel) || WithinAbs(ref_a, tol_abs));
    CHECK_THAT(compute_molar_heat_capacity_v(T_b, debye_b, n_b),
      WithinRel(ref_b, tol_rel) || WithinAbs(ref_b, tol_abs));
  }
  SECTION("compute_helmholtz_free_energy") {
    double ref_a = -173316.33430428786;
    double ref_b = -1445.5499830398378;
    CHECK_THAT(compute_helmholtz_free_energy(T_a, debye_a, n_a),
      WithinRel(ref_a, tol_rel) || WithinAbs(ref_a, tol_abs));
    CHECK_THAT(compute_helmholtz_free_energy(T_b, debye_b, n_b),
      WithinRel(ref_b, tol_rel) || WithinAbs(ref_b, tol_abs));
  }
  SECTION("compute_entropy") {
    double ref_a = 131.64954770825574;
    double ref_b = 13.356996569280064;
    CHECK_THAT(compute_entropy(T_a, debye_a, n_a),
      WithinRel(ref_a, tol_rel) || WithinAbs(ref_a, tol_abs));
    CHECK_THAT(compute_entropy(T_b, debye_b, n_b),
      WithinRel(ref_b, tol_rel) || WithinAbs(ref_b, tol_abs));
  }
  SECTION("compute_dmolar_heat_capacity_v_dT") {
    double ref_a = 0.0001828985485514273;
    double ref_b = 0.035412634394977965;
    CHECK_THAT(compute_dmolar_heat_capacity_v_dT(T_a, debye_a, n_a),
      WithinRel(ref_a, tol_rel) || WithinAbs(ref_a, tol_abs));
    CHECK_THAT(compute_dmolar_heat_capacity_v_dT(T_b, debye_b, n_b),
      WithinRel(ref_b, tol_rel) || WithinAbs(ref_b, tol_abs));
  }
}
