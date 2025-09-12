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
#include "burnman/optim/roots/bracket.hpp"
#include "burnman/optim/roots/brent.hpp"
#include "tolerances.hpp"

using namespace Catch::Matchers;
using namespace burnman;
using namespace optim::roots;

// Simple test function
double test_function(double x, void*) {
    return x * x - 4; // Roots at x = -2 and x = 2
}

// Simple test function with parameter(s)
struct test_params {
  double y;
};
double test_function_with_p(double x, void* p) {
  auto* params = static_cast<const test_params*>(p);
  return x * x - (4.0 + 10.0 * params->y);
}

TEST_CASE("Test brent method", "[optim][roots]") {
  SECTION("Root at x = 2") {
    double x_lo = 1.0;
    double x_hi = 3.0;
    double root;
    REQUIRE_NOTHROW(root = brent(&test_function, nullptr, x_lo, x_hi));
    REQUIRE_THAT(root, WithinRel(2.0, tol_rel) || WithinAbs(2.0, tol_abs));
  }
  SECTION("Root at x = -2") {
    double x_lo = -3.0;
    double x_hi = -1.0;
    double root;
    REQUIRE_NOTHROW(root = brent(&test_function, nullptr, x_lo, x_hi));
    REQUIRE_THAT(root, WithinRel(-2.0, tol_rel) || WithinAbs(-2.0, tol_abs));
  }
  SECTION("No root") {
    double x_lo = 5.0;
    double x_hi = 6.0;
    double root;
    REQUIRE_THROWS(root = brent(&test_function, nullptr, x_lo, x_hi));
  }
  SECTION("Root at edge") {
    double x_lo = 2.0;
    double x_hi = 3.0;
    double root;
    REQUIRE_NOTHROW(root = brent(&test_function, nullptr, x_lo, x_hi));
    REQUIRE_THAT(root, WithinRel(2.0, tol_rel) || WithinAbs(2.0, tol_abs));
  }
  SECTION("Function with parameter") {
    double x_lo = 2.5;
    double x_hi = 3.5;
    test_params p{0.5};
    double root;
    REQUIRE_NOTHROW(root = brent(&test_function_with_p, p, x_lo, x_hi));
    REQUIRE_THAT(root, WithinRel(3.0, tol_rel) || WithinAbs(3.0, tol_abs));
  }
}

TEST_CASE("Test bracket_root", "[optim][roots]") {
  SECTION("Test with valid bracket") {
    double x_lo = 1.0;
    double x_hi = 3.0;
    REQUIRE(bracket_root(&test_function, nullptr, x_lo, x_hi, 0.1));
  }
  SECTION("Test bracket expansion - downhill low") {
    double x_lo = 3.0;
    double x_hi = 3.5;
    REQUIRE(bracket_root(&test_function, nullptr, x_lo, x_hi, 0.1));
    REQUIRE(x_lo < 3.0);
    REQUIRE(x_hi == 3.5);
  }
  SECTION("Test bracket expansion - downhill high") {
    double x_lo = 1.0;
    double x_hi = 1.5;
    REQUIRE(bracket_root(&test_function, nullptr, x_lo, x_hi, 0.1));
    REQUIRE(x_lo == 1.0);
    REQUIRE(x_hi > 1.5);
  }
  SECTION("Test invalid bracket") {
    double x_lo = 10.0;
    double x_hi = 11.0;
    REQUIRE_FALSE(bracket_root(&test_function, nullptr, x_lo, x_hi, 0.1, 1.0, 10));
  }
  SECTION("Test root at edge") {
    double x_lo = 2.0;
    double x_hi = 3.0;
    REQUIRE(bracket_root(&test_function, nullptr, x_lo, x_hi, 0.1));
  }
  SECTION("Test expansion factor > 1") {
    double x_lo = 3.0;
    double x_hi = 3.5;
    // With expansion factor of 10, should bracket root in 2 iterations
    REQUIRE(bracket_root(&test_function, nullptr, x_lo, x_hi, 0.1, 10.0, 2));
  }
  SECTION("Test function with parameters") {
    double x_lo = 3.0;
    double x_hi = 3.5;
    SECTION("With root in interval") {
      test_params p_with_root{0.0};
      REQUIRE(bracket_root(&test_function_with_p, p_with_root, x_lo, x_hi, 0.1));  
    }
    SECTION("With no root in interval") {
      test_params p_no_root{20.0}; // No root in 100 iterations with dx = 0.1
      REQUIRE_FALSE(bracket_root(&test_function_with_p, p_no_root, x_lo, x_hi, 0.1));
    }
  }
}
