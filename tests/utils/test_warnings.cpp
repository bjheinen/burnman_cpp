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
#include "burnman/utils/warnings.hpp"
#include <iostream>
#include <sstream>
#include <string>

using namespace burnman;

TEST_CASE("Validate warnings", "[utils][warnings]") {
  // Temporarily allow warnings
  utils::suppress_warnings = false;
  // Redirect cerr to grab output
  std::streambuf* orig_cerr = std::cerr.rdbuf();
  std::ostringstream oss;
  std::cerr.rdbuf(oss.rdbuf());
  // Generate a warning
  utils::warn("This is a test warning!");
  REQUIRE(oss.str() == "Warning: This is a test warning!\n");
  // Restore cerr buffer
  std::cerr.rdbuf(orig_cerr);
  // Block warnings again
  utils::suppress_warnings = true;
}
