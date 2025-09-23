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
#include "burnman/utils/validate_optionals.hpp"
#include <optional>
#include <stdexcept>
#include <sstream>
#include <string>
#include "burnman/utils/warnings.hpp"
#include "tolerances.hpp"

using namespace Catch::Matchers;
using namespace burnman;

TEST_CASE("require_set", "[utils][validate_optionals]") {
  SECTION("T double") {
    std::optional<double> opt;
    REQUIRE_THROWS_AS(utils::require_set(opt, "Test parameter"), std::invalid_argument);
    opt = 42.0;
    REQUIRE_NOTHROW(utils::require_set(opt, "Test parameter"));
  }
  SECTION("T string") {
    std::optional<std::string> opt;
    REQUIRE_THROWS_AS(utils::require_set(opt, "Test parameter"), std::invalid_argument);
    opt = "Test";
    REQUIRE_NOTHROW(utils::require_set(opt, "Test parameter"));
  }
}

TEST_CASE("fallback_to_default", "[utils][validate_optionals]") {
  SECTION("T double") {
    std::optional<double> opt;
    utils::fallback_to_default(opt, 3.14, false, "Test parameter");
    REQUIRE(opt.has_value());
    REQUIRE_THAT(*opt, WithinRel(3.14, tol_rel) || WithinAbs(3.14, tol_abs));
    utils::fallback_to_default(opt, 2.71, false, "Test parameter");
    REQUIRE(opt.has_value());
    REQUIRE_THAT(*opt, WithinRel(3.14, tol_rel) || WithinAbs(3.14, tol_abs));
  }
  SECTION("T string") {
    std::optional<std::string> opt;
    utils::fallback_to_default(opt, "default", false, "Test parameter");
    REQUIRE(opt.has_value());
    REQUIRE(*opt == "default");
    utils::fallback_to_default(opt, "new default", false, "Test parameter");
    REQUIRE(opt.has_value());
    REQUIRE(*opt == "default");
  }
}

TEST_CASE("check_in_range", "[utils][validate_optionals]") {
  // Suppress warnings for test
  std::optional<double> opt;
  opt = 5.0;
  REQUIRE(utils::check_in_range(opt, 0.0, 10.0));
  REQUIRE_FALSE(utils::check_in_range(opt, 6.0, 10.0));
  REQUIRE_THROWS_AS(utils::check_in_range(opt, 6.0, 10.0, "Test parameter", true), std::runtime_error);
}
