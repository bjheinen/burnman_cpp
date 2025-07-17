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
#include "burnman/utils/constants.hpp"

// Compile time checks that physical constants haven't been modified
TEST_CASE("constants::physics CODATA 2022 values", "[core][utils][constants]") {
  STATIC_REQUIRE(constants::physics::gas_constant == 8.314462618);
  STATIC_REQUIRE(constants::physics::boltzmann == 1.380649e-23);
  STATIC_REQUIRE(constants::physics::avogadro == 6.02214076e23);
  STATIC_REQUIRE(constants::physics::dirac == 1.054571817e-34);
  STATIC_REQUIRE(constants::physics::gravitation == 6.67430e-11);
}
