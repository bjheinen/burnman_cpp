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
#include "burnman/utils/math_utils.hpp"
#include <cmath>
#include <Eigen/Dense>
#include "tolerances.hpp"

using namespace Catch::Matchers;
using namespace burnman;

TEST_CASE("logish is log for values > eps; 1D", "[utils][math_utils]") {
  Eigen::ArrayXd x(3);
  x << 1.0e-6, 1.0e-5, 1.0e-3; // All > eps = 1e-7
  Eigen::ArrayXd expected = x.log();
  Eigen::ArrayXd result = utils::logish(x);
  CHECK(result.isApprox(expected, tol_rel));
}

TEST_CASE("logish is log for values > eps; 2D", "[utils][math_utils]") {
  Eigen::ArrayXXd x(2, 2);
  x <<
    1.0e-6, 1.0e-5,
    2.0e1, 1.0e-3; // All > eps = 1e-7
  Eigen::ArrayXXd expected = x.log();
  Eigen::ArrayXXd result = utils::logish(x);
  CHECK(result.isApprox(expected, tol_rel));
}

TEST_CASE("logish uses Taylor expansion below eps", "[utils][math_utils]") {
  Eigen::ArrayXd x(6);
  x << 2.0e-8, 1.4e-9, 9.999e-8, 5.0e-13, 1.0e-35, 0;  // All <= eps (1.e-7)
  Eigen::ArrayXd result = utils::logish(x);
  double lim = std::log(1.0e-7) - 1.5;
  // log(x) ~ logish(x) at x ~ eps
  double ref_a = std::log(x(2));
  CHECK_THAT(result(2),
    WithinRel(ref_a, tol_rel) || WithinAbs(ref_a, tol_abs));
  // logish(x) same order of magnitude at near eps
  CHECK(int(result(0)) == int(std::log(x(0))));
  // logish levels
  REQUIRE(int(std::log(x(1))) == -20);
  CHECK(int(result(1)) == -17);
  REQUIRE(int(std::log(x(3))) == -28);
  CHECK(int(result(3)) == -17);
  // log(x->0) -> inf; logish(x->0) -> lim (log(1e-7)-1.5)
  CHECK_THAT(result(4),
    WithinRel(lim, tol_rel) || WithinAbs(lim, tol_abs));
  CHECK_THAT(result(5),
    WithinRel(lim, tol_rel) || WithinAbs(lim, tol_abs));
}

TEST_CASE("logish - mixed range 2D", "[utils][math_utils]") {
  Eigen::ArrayXXd x(2, 3);
  x <<
    2.0e-1, 1.4e-5, 9.99e-8,
    5.0e-13, 1.0e1, 0;
  Eigen::ArrayXXd result;
  REQUIRE_NOTHROW(result = utils::logish(x));
  REQUIRE(std::isfinite(result(1,2)));
  REQUIRE(result.rows() == 2);
  REQUIRE(result.cols() == 3);
}

TEST_CASE("inverseish is 1/x for values > eps; 1D", "[utils][math_utils]") {
  Eigen::ArrayXd x(3);
  x << 1.0e-4, 1.0e-2, 1.0e-4; // All > eps = 1e-5
  Eigen::ArrayXd expected = 1.0 / x;
  Eigen::ArrayXd result = utils::inverseish(x);
  CHECK(result.isApprox(expected, tol_rel));
}

TEST_CASE("inverseish is 1/x for values > eps; 2D", "[utils][math_utils]") {
  Eigen::ArrayXXd x(2, 2);
  x <<
    1.0e-4, 1.0e-2,
    2.0e1, 1.0e-3; // All > eps = 1e-5
  Eigen::ArrayXXd expected = x.inverse();
  Eigen::ArrayXXd result = utils::inverseish(x);
  CHECK(result.isApprox(expected, tol_rel));
}

TEST_CASE("inverseish uses Taylor expansion below eps", "[utils][math_utils]") {
  Eigen::ArrayXd x(6);
  x << 2.0e-6, 1.4e-9, 9.999e-6, 5.0e-13, 1.0e-35, 0;  // All <= eps
  Eigen::ArrayXd result = utils::inverseish(x);
  Eigen::ArrayXd ref(6);
  ref << 180000, 199986, 100010, 199999.995, 200000, 200000;
  REQUIRE(result.isApprox(ref, tol_rel));
}

TEST_CASE("inverseish - mixed range 2D", "[utils][math_utils]") {
  Eigen::ArrayXXd x(2, 3);
  x <<
    2.0e-1, 1.4e-5, 9.99e-8,
    5.0e-13, 1.0e1, 0;
  Eigen::ArrayXXd result;
  REQUIRE_NOTHROW(result = utils::inverseish(x));
  REQUIRE(std::isfinite(result(1,2)));
  REQUIRE(result.rows() == 2);
  REQUIRE(result.cols() == 3);
}
