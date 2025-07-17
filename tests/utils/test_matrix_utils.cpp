/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "tolerances.hpp"
#include "burnman/utils/matrix_utils.hpp"

TEST_CASE("jagged2square n=5", "[utils][matrix_utils]") {
  int n = 5;
  std::vector<std::vector<double>> v = {
    {0.0, 24.74e3, 26.0e3, 24.3e3},
    {24.74e3, 0.0, 0.0e3},
    {60.53136e3, 0.0},
    {10.0e3}
  };  
  Eigen::MatrixXd expected(n, n);
  expected <<
    0, 0, 24740,    26000, 24300,
    0, 0, 24740,        0,     0,
    0, 0,     0, 60531.36,     0,
    0, 0,     0,        0, 10000,
    0, 0,     0,        0,     0;
  Eigen::MatrixXd result = utils::jagged2square(v, n);
  REQUIRE((result.array() == expected.array()).all());
  REQUIRE(result.isUpperTriangular());
}

TEST_CASE("jagged2square n=5 (equal rows)", "[utils][matrix_utils]") {
  int n = 5;
  std::vector<std::vector<double>> v = {
    {0.0, 24.74e3, 26.0e3, 24.3e3},
    {24.74e3, 0.0, 0.0e3},
    {60.53136e3, 0.0},
    {10.0e3},
    {}
  };  
  Eigen::MatrixXd expected(n, n);
  expected <<
    0, 0, 24740,    26000, 24300,
    0, 0, 24740,        0,     0,
    0, 0,     0, 60531.36,     0,
    0, 0,     0,        0, 10000,
    0, 0,     0,        0,     0;
  Eigen::MatrixXd result = utils::jagged2square(v, n);
  REQUIRE((result.array() == expected.array()).all());
  REQUIRE(result.isUpperTriangular());
}

TEST_CASE("jagged2square n=2", "[utils][matrix_utils]") {
  int n = 2;
  std::vector<std::vector<double>> v = {
    {18.0e3}
  };
  Eigen::MatrixXd expected(n, n);
  expected <<
    0, 18.0e3,
    0,      0;
  Eigen::MatrixXd result = utils::jagged2square(v, n);
  REQUIRE((result.array() == expected.array()).all());
  REQUIRE(result.isUpperTriangular());
}

TEST_CASE("jagged2square (empty input)") {
  std::vector<std::vector<double>> v;
  Eigen::MatrixXd result = utils::jagged2square(v, 3);
  REQUIRE(result.isZero());
  REQUIRE(result.isUpperTriangular());
}
