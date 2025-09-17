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
#include "burnman/utils/types/ndarray.hpp"
#include <cstddef>
#include <vector>
#include "burnman/utils/vector_utils.hpp"

using namespace burnman;

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

TEST_CASE("NDArray construction", "[utils][types]") {
  SECTION("Default constructor") {
    types::NDArray<int> arr;
    REQUIRE(arr.shape().empty());
    REQUIRE(arr.size() == 0);
    REQUIRE(arr.data().empty());
  }
  SECTION("Construction with shape") {
    types::NDArray<double> arr({2, 3});
    REQUIRE(arr.shape() == std::vector<std::size_t>{2, 3});
    REQUIRE(arr.size() == 6);
    REQUIRE(arr.data().size() == 6);
    for (auto v : arr.data()) {
      REQUIRE(v == double{});
    }
  }
}

TEST_CASE("NDArray indexing", "[utils][types]") {
  types::NDArray<int> arr({2, 4});
  REQUIRE(arr.size() == 8);
  SECTION("Flat indexing - [] & ()") {
    REQUIRE_NOTHROW(arr[1] = 10);
    REQUIRE_NOTHROW(arr(2) = 20);
    REQUIRE(arr(1) == 10);
    REQUIRE(arr[2] == 20);
    const types::NDArray<int>& const_arr = arr;
    REQUIRE(const_arr(1) == 10);
    REQUIRE(const_arr[2] == 20);
  }
  SECTION("Grid indexing - vector & init_list") {
    std::vector<std::size_t> ind_0{0, 0};
    std::vector<std::size_t> ind_n{1, 2};
    REQUIRE_NOTHROW(arr({0, 0}) = 1);
    REQUIRE_NOTHROW(arr(ind_n) = 2);
    REQUIRE(arr(ind_0) == 1);
    REQUIRE(arr({1, 2}) == 2);
    const types::NDArray<int>& const_arr = arr;
    REQUIRE(const_arr({ind_0}) == 1);
    REQUIRE(const_arr({ind_n}) == 2);
  }
  SECTION("set_shape") {
    arr({0, 0}) = 42;
    REQUIRE(arr.size() == 8);
    arr.set_shape({3, 3});
    REQUIRE(arr.shape() == std::vector<std::size_t>{3, 3});
    REQUIRE(arr.size() == 9);
    REQUIRE(arr.data().size() == 9);
    REQUIRE_FALSE(arr[0] == 42);
    REQUIRE(arr[0] == int{});
  }
}
TEST_CASE("NDArray vs utils::vector_utils", "[utils][types]") {
  std::vector<std::size_t> shape{2, 3, 4};
  std::vector<std::size_t> strides = utils::compute_strides(shape);
  types::NDArray<int> arr(shape);
  for (std::size_t i = 0; i < arr.size(); ++ i) {
    arr[i] = static_cast<int>(i);
  }
  for (std::size_t i = 0; i < arr.size(); ++i) {
    auto indices = utils::map_index(i, strides);
    REQUIRE(arr(indices) == static_cast<int>(i));
  }
}
