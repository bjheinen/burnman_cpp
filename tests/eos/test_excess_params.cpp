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
#include "burnman/eos/components/excess_params.hpp"

using namespace burnman;

TEST_CASE("Test excesses::Excesses", "[utils][types]") {
  eos::excesses::Excesses e;
  SECTION("Default values") {
    REQUIRE(e.G == 0.0);
    REQUIRE(e.dGdT == 0.0);
    REQUIRE(e.dGdP == 0.0);
    REQUIRE(e.d2GdT2 == 0.0);
    REQUIRE(e.d2GdP2 == 0.0);
    REQUIRE(e.d2GdPdT == 0.0);
  }
  SECTION("operator +=") {
    e.G = 1.0;
    e.dGdT = 2.0;
    e.dGdP = 3.0;
    eos::excesses::Excesses e2;
    e2.G = 0.5;
    e2.dGdT = -1.0;
    e2.dGdP = 4.0;
    e2.d2GdT2 = 10.0;
    e2.d2GdP2 = 20.0;
    e2.d2GdPdT = 30.0;
    e += e2;
    REQUIRE(e.G == 1.5);
    REQUIRE(e.dGdT == 1.0);
    REQUIRE(e.dGdP == 7.0);
    REQUIRE(e.d2GdT2 == 10.0);
    REQUIRE(e.d2GdP2 == 20.0);
    REQUIRE(e.d2GdPdT == 30.0);
  }
}
