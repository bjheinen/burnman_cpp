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
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "burnman/tools/averaging/averaging_schemes.hpp"
#include "tolerances.hpp"

using namespace Catch::Matchers;
using namespace burnman;
using namespace averaging;

struct AveragingSchemeFixture {
  Eigen::ArrayXd v_a = (Eigen::ArrayXd(2) << 0.25, 0.75).finished();
  Eigen::ArrayXd X_a = (Eigen::ArrayXd(2) << 1.0, 2.0).finished();
  Eigen::ArrayXd v_b = (Eigen::ArrayXd(4) << 0.3, 0.45, 0.2, 0.05).finished();
  Eigen::ArrayXd X_b = (Eigen::ArrayXd(4) << 120.0e9, 208.45e9, 264.2e9, 76.1e9).finished();
  Eigen::ArrayXd X_c = (Eigen::ArrayXd(4) << 96.0e9, 173.2e9, 207.35, 36.2).finished();
  double ref_voigt = 186.4475e9;
  double ref_reuss = 164668047448.38943;
  double ref_vrh = 175557773724.1947;
  double ref_hsl_K = 164668047455.9929;
  double ref_hsl_G = 705.0786544800931;
  double ref_hsu_K = 177929798385.34683;
  double ref_hsu_G = 86493747346.95517;
};

template <typename T>
struct AveragingSchemeMulti : AveragingSchemeFixture {
  T scheme;
};

TEST_CASE_METHOD(AveragingSchemeFixture, "Averaging::Voigt", "[core][averaging_schemes]") {
  // Quick check of Voigt average
  REQUIRE(utils::voigt_fn(v_a, X_a) == 1.75);
  // Check voigt_fn normalises if needed
  REQUIRE(utils::voigt_fn(v_a*34.235, X_a) == 1.75);
  CHECK_THAT(utils::voigt_fn(v_b, X_b),
    WithinRel(ref_voigt, tol_rel) ||
    WithinAbs(ref_voigt, tol_abs));
  // Compare to scheme interface
  Voigt scheme;
  REQUIRE(scheme.average_bulk_moduli(v_b, X_b, X_c)
    == scheme.average_shear_moduli(v_b, X_c, X_b));
  REQUIRE(scheme.average_bulk_moduli(v_b, X_b, X_c) == utils::voigt_fn(v_b, X_b));
}

TEST_CASE_METHOD(AveragingSchemeFixture, "Averaging::Reuss", "[core][averaging_schemes]") {
  // Quick check of Reuss average
  REQUIRE(utils::reuss_fn(v_a, X_a) == 1.6);
  // Check reuss_fn normalises if needed
  REQUIRE(utils::reuss_fn(v_a*34.235, X_a) == 1.6);
  CHECK_THAT(utils::reuss_fn(v_b, X_b),
    WithinRel(ref_reuss, tol_rel) ||
    WithinAbs(ref_reuss, tol_abs));
  // Compare to scheme interface
  Reuss scheme;
  REQUIRE(scheme.average_bulk_moduli(v_b, X_b, X_c)
    == scheme.average_shear_moduli(v_b, X_c, X_b));
  REQUIRE(scheme.average_bulk_moduli(v_b, X_b, X_c) == utils::reuss_fn(v_b, X_b));
}

TEST_CASE_METHOD(AveragingSchemeFixture, "Averaging::VoigtReussHill", "[core][averaging_schemes]") {
  // Quick check of VRH average
  REQUIRE(utils::voigt_reuss_hill_fn(v_a, X_a) == 1.675);
  CHECK_THAT(utils::voigt_reuss_hill_fn(v_b, X_b),
    WithinRel(ref_vrh, tol_rel) ||
    WithinAbs(ref_vrh, tol_abs));
  // Compare to scheme interface
  VoigtReussHill scheme;
  REQUIRE(scheme.average_bulk_moduli(v_b, X_b, X_c)
    == scheme.average_shear_moduli(v_b, X_c, X_b));
  REQUIRE(scheme.average_bulk_moduli(v_b, X_b, X_c) == utils::voigt_reuss_hill_fn(v_b, X_b));
}

TEMPLATE_TEST_CASE_METHOD(AveragingSchemeMulti, "Averaging (common fns)",
  "[core][averaging_schemes]",
  Voigt, Reuss, VoigtReussHill,
  HashinShtrikmanLower, HashinShtrikmanUpper, HashinShtrikman
) {
  auto& sch = this->scheme;
  auto& va = this->v_a;
  auto& Xa = this->X_a;
  CHECK(sch.average_density(va, Xa) == 1.75);
  CHECK(sch.average_density(va*2.0, Xa) == 1.75);
  CHECK(sch.average_thermal_expansivity(va, Xa) == 1.75);
  CHECK(sch.average_thermal_expansivity(va*2.0, Xa) == 1.75);
  CHECK(sch.average_heat_capacity_v(va, Xa) == 1.75);
  CHECK(sch.average_heat_capacity_p(va, Xa) == 1.75);
}

TEST_CASE_METHOD(AveragingSchemeFixture, "Averaging::HashinShtrikmanLower", "[core][averaging_schemes]") {
  HashinShtrikmanLower scheme;
  CHECK_THAT(scheme.average_bulk_moduli(v_b, X_b, X_c),
    WithinRel(ref_hsl_K, tol_rel) ||
    WithinAbs(ref_hsl_K, tol_abs));
  CHECK_THAT(scheme.average_shear_moduli(v_b, X_b, X_c),
    WithinRel(ref_hsl_G, tol_rel) ||
    WithinAbs(ref_hsl_G, tol_abs));
}

TEST_CASE_METHOD(AveragingSchemeFixture, "Averaging::HashinShtrikmanUpper", "[core][averaging_schemes]") {
  HashinShtrikmanUpper scheme;
  CHECK_THAT(scheme.average_bulk_moduli(v_b, X_b, X_c),
    WithinRel(ref_hsu_K, tol_rel) ||
    WithinAbs(ref_hsu_K, tol_abs));
  CHECK_THAT(scheme.average_shear_moduli(v_b, X_b, X_c),
    WithinRel(ref_hsu_G, tol_rel) ||
    WithinAbs(ref_hsu_G, tol_abs));
}

TEST_CASE_METHOD(AveragingSchemeFixture, "Averaging::HashinShtrikman", "[core][averaging_schemes]") {
  HashinShtrikman scheme;
  double ref_hs_K = (ref_hsl_K + ref_hsu_K) / 2.0;
  double ref_hs_G = (ref_hsl_G + ref_hsu_G) / 2.0;
  CHECK_THAT(scheme.average_bulk_moduli(v_b, X_b, X_c),
    WithinRel(ref_hs_K, tol_rel) ||
    WithinAbs(ref_hs_K, tol_abs));
  CHECK_THAT(scheme.average_shear_moduli(v_b, X_b, X_c),
    WithinRel(ref_hs_G, tol_rel) ||
    WithinAbs(ref_hs_G, tol_abs));
}
