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
#include "burnman/utils/vector_utils.hpp"
#include <cstddef>
#include <vector>
#include "tolerances.hpp"

using namespace burnman;

TEST_CASE("compute_strides", "[utils][vector_utils]") {
  SECTION("1D shape") {
    std::vector<std::size_t> shape{5};
    std::vector<std::size_t> strides = utils::compute_strides(shape);
    REQUIRE(strides == std::vector<std::size_t>{1});
  }
  SECTION("2D shape") {
    std::vector<std::size_t> shape{2, 4};
    std::vector<std::size_t> strides = utils::compute_strides(shape);
    REQUIRE(strides == std::vector<std::size_t>{4, 1});
  }
  SECTION("3D shape") {
    std::vector<std::size_t> shape{3, 4, 5};
    std::vector<std::size_t> strides = utils::compute_strides(shape);
    REQUIRE(strides == std::vector<std::size_t>{20, 5, 1});
  }
  SECTION("empty") {
    std::vector<std::size_t> shape{};
    std::vector<std::size_t> strides = utils::compute_strides(shape);
    REQUIRE(strides.empty());
  }
}

TEST_CASE("flatten_index", "[utils][vector_utils]") {
  SECTION("2D shape") {
    std::vector<std::size_t> shape{2, 4};
    std::vector<std::size_t> strides = utils::compute_strides(shape);
    REQUIRE(utils::flatten_index({0, 0}, strides) == 0);
    REQUIRE(utils::flatten_index({1, 3}, strides) == 7);
  }
  SECTION("3D shape") {
    std::vector<std::size_t> shape{3, 4, 5};
    std::vector<std::size_t> strides = utils::compute_strides(shape);
    REQUIRE(utils::flatten_index({0, 0, 0}, strides) == 0);
    REQUIRE(utils::flatten_index({0, 1, 0}, strides) == 5);
    REQUIRE(utils::flatten_index({2, 3, 4}, strides) == 59);
  }
  SECTION("bad input") {
    std::vector<std::size_t> strides{4, 1};
    REQUIRE_THROWS_AS(utils::flatten_index({0}, strides), std::invalid_argument);
    REQUIRE_THROWS_AS(utils::flatten_index({0, 1, 2}, strides), std::invalid_argument);
  }
}

TEST_CASE("map_index", "[utils][vector_utils]") {
  SECTION("2D shape") {
    std::vector<std::size_t> shape{2, 4};
    std::vector<std::size_t> strides = utils::compute_strides(shape);
    REQUIRE(utils::map_index(0, strides) == std::vector<std::size_t>{0, 0});
    REQUIRE(utils::map_index(7, strides) == std::vector<std::size_t>{1, 3});
  }
  SECTION("3D shape") {
    std::vector<std::size_t> shape{3, 4, 5};
    std::vector<std::size_t> strides = utils::compute_strides(shape);
    REQUIRE(utils::map_index(0, strides) == std::vector<std::size_t>{0, 0, 0});
    REQUIRE(utils::map_index(5, strides) == std::vector<std::size_t>{0, 1, 0});
    REQUIRE(utils::map_index(59, strides) == std::vector<std::size_t>{2, 3, 4});
  }
}

TEST_CASE("map_index / flatten_index check", "[utils][vector_utils]") {
  std::vector<std::size_t> shape{3, 5, 6, 2};
  std::vector<std::size_t> strides = utils::compute_strides(shape);
  for (std::size_t flat = 0; flat < 180; flat++) {
    std::vector<std::size_t> indices = utils::map_index(flat, strides);
    std::size_t flat_calc = utils::flatten_index(indices, strides);
    REQUIRE(flat == flat_calc);
  }
}
