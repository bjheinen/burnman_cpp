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
#include "burnman/utils/string_utils.hpp"
#include <string>
#include <vector>
#include "tolerances.hpp"

TEST_CASE("Extract numeric prefix", "[utils][string_utils]") {
  REQUIRE(utils::extract_numeric_prefix("2SiO2") == "2");
  REQUIRE(utils::extract_numeric_prefix("0.45SiO2") == "0.45");
  REQUIRE(utils::extract_numeric_prefix("1/2Fe2O3") == "1/2");
  REQUIRE(utils::extract_numeric_prefix("Al2O3") == "");
  REQUIRE(utils::extract_numeric_prefix("2e1Fe") == "2e1");
  REQUIRE(utils::extract_numeric_prefix("2e-1Fe") == "2e-1");
  // Doesn't work with E!
  REQUIRE(utils::extract_numeric_prefix("2E-1Fe") != "2E-1");
}

TEST_CASE("Extended stod", "[utils][string_utils]") {
  REQUIRE(utils::stod("2") == 2);
  REQUIRE(utils::stod("0.45") == 0.45);
  REQUIRE(utils::stod("1/2") == 0.5);
  REQUIRE(utils::stod("1e2") == 100);
}

TEST_CASE("Utils string join", "[utils][string_utils]") {
  std::vector<std::string> s1 = {"A", "B"};
  std::vector<std::string> s2 = {"A", "B", "C"};
  REQUIRE(utils::join(s1) == "AB");
  REQUIRE(utils::join(s2) == "ABC");
  REQUIRE(utils::join(s1, " ") == "A B");
  REQUIRE(utils::join(s2, ", ") == "A, B, C");
}
