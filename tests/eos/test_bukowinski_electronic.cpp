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
#include "burnman/eos/bukowinski_electronic.hpp"
#include "tolerances.hpp"

using namespace bukowinski;
using namespace Catch::Matchers;

TEST_CASE("Check zero returns", "[eos][bukowinski]") {
  MineralParams params;
  params.V_0 = 7.0e-6;
  params.T_0 = 300.0;
  params.bel_0 = 0.004;
  params.gel = 1.5;
  SECTION("Reference T returns 0") {
    double T = 300.0;
    double V = 6.0e-6;
    REQUIRE(compute_helmholtz_el(T, V, params) == 0);
    REQUIRE(compute_pressure_el(T, V, params) == 0);
    REQUIRE(compute_KT_over_V(T, V, params) == 0);
  }
}

// Reference values from Py burnman v2.1.1a0
TEST_CASE("Bukowinski functions python reference values", "[eos][bukowinski]") {
  MineralParams params;
  params.V_0 = 7.0e-6;
  params.T_0 = 300.0;
  params.bel_0 = 0.004;
  params.gel = 1.5;
  double T = 1200.0;
  double V = 2e-6;

  SECTION("compute_helmholtz_el") {
    double ref = -412.3459160934548;
    CHECK_THAT(compute_helmholtz_el(T, V, params),
      WithinRel(ref, tol_rel) || WithinAbs(ref, tol_abs));
  }
  SECTION("compute_pressure_el") {
    double ref = 309259437.0700911;
    CHECK_THAT(compute_pressure_el(T, V, params),
      WithinRel(ref, tol_rel) || WithinAbs(ref, tol_abs));
  }
  SECTION("compute_entropy_el") {
    double ref = 0.733059406388364;
    CHECK_THAT(compute_entropy_el(T, V, params),
      WithinRel(ref, tol_rel) || WithinAbs(ref, tol_abs));
  }
  SECTION("compute_KT_over_V") {
    double ref = -77314859267522.78;
    CHECK_THAT(compute_KT_over_V(T, V, params),
      WithinRel(ref, tol_rel) || WithinAbs(ref, tol_abs));
  }
  SECTION("compute_CV_over_T") {
    double ref = 0.00061088283865697;
    CHECK_THAT(compute_CV_over_T(V, params),
      WithinRel(ref, tol_rel) || WithinAbs(ref, tol_abs));
  }
  SECTION("compute_alpha_KT") {
    double ref = 549794.5547912731;
    CHECK_THAT(compute_alpha_KT(T, V, params),
      WithinRel(ref, tol_rel) || WithinAbs(ref, tol_abs));
  }
}
