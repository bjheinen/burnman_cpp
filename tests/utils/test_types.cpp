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
#include "burnman/utils/types/simple_types.hpp"

//TODO: Update directory structure - test new types

TEST_CASE("FormulaMap addition and scalar multiplication", "[utils][types]") {
  
  types::FormulaMap fm1 {{"Si", 1.0}, {"O", 2.0}};
  types::FormulaMap fm2 {{"Mg", 1.0}, {"O", 1.0}};

  SECTION("Addition") {
    types::FormulaMap result = fm1 + fm2;
    REQUIRE(result["Si"] == 1.0);
    REQUIRE(result["Mg"] == 1.0);
    REQUIRE(result["O"] == 3.0);
  }

  SECTION("In-place addition") {
    fm1 += fm2;
    REQUIRE(fm1["Si"] == 1.0);
    REQUIRE(fm1["Mg"] == 1.0);
    REQUIRE(fm1["O"] == 3.0);
  }

  SECTION("Scalar multiplication") {
    types::FormulaMap result = fm1 * 2.0;
    REQUIRE(result["Si"] == 2.0);
    REQUIRE(result["O"] == 4.0);
  }

  SECTION("Scalar multiplication (swapped order)") {
    types::FormulaMap result = 3.0 * fm1;
    REQUIRE(result["Si"] == 3.0);
    REQUIRE(result["O"] == 6.0);
  }

  SECTION("In-place scalar multiplication") {
    fm1 *= 2.5;
    REQUIRE(fm1["Si"] == 2.5);
    REQUIRE(fm1["O"] == 5.0);
  }
}

// TODO: test excesses::Excesses (and move to types.hpp)
