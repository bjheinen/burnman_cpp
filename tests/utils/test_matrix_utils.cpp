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
  Eigen::Index n = 5;
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
  Eigen::Index n = 5;
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
  Eigen::Index n = 2;
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

TEST_CASE("jagged2square (empty input)", "[utils][matrix_utils]") {
  std::vector<std::vector<double>> v;
  Eigen::MatrixXd result = utils::jagged2square(v, 3);
  REQUIRE(result.isZero());
  REQUIRE(result.isUpperTriangular());
}

TEST_CASE("populate_interaction_matrix; check shape", "[utils][matrix_utils]") {
  Eigen::Index n = 3;
  Eigen::ArrayXd alphas(n);
  alphas << 1.0, 2.0, 3.0;
  Eigen::MatrixXd interaction(n, n);
  interaction << 1, 2, 3,
                 4, 5, 6,
                 7, 8, 9;
  Eigen::MatrixXd result = utils::populate_interaction_matrix(interaction, alphas, n);
  REQUIRE(result.rows() == n);
  REQUIRE(result.cols() == n);
  REQUIRE(result.isUpperTriangular());
}

TEST_CASE("populate_interaction_matrix; values", "[utils][matrix_utils]") {
  Eigen::Index n = 5;
  Eigen::MatrixXd interaction(n, n);
  interaction <<
    0, 0, 24740,    26000, 24300,
    0, 0, 24740,        0,     0,
    0, 0,     0, 60531.36,     0,
    0, 0,     0,        0, 10000,
    0, 0,     0,        0,     0;

  SECTION("Symmetric ((alphas == 1).all())") {
    Eigen::ArrayXd alphas(n);
    alphas << 1.0, 1.0, 1.0, 1.0, 1.0;
    Eigen::MatrixXd result = utils::populate_interaction_matrix(interaction, alphas, n);
    REQUIRE((result.array() == interaction.array()).all());
  }
  SECTION("Asymmetric (!alphas == 1).all())") {
    Eigen::ArrayXd alphas(n);
    alphas << 1, 2, 3, 4, 5;
    Eigen::MatrixXd expected(n, n);
    expected <<
      0, 0, 24740*0.5,          26000*0.4, 24300*(1.0/3.0),
      0, 0, 24740*0.4,                  0,               0,
      0, 0,         0, 60531.36*(2.0/7.0),               0,
      0, 0,         0,                  0, 10000*(2.0/9.0),
      0, 0,         0,                  0,               0;
    Eigen::MatrixXd result = utils::populate_interaction_matrix(interaction, alphas, n);
    REQUIRE(result.isApprox(expected, tol_rel));
  }
  SECTION("Zero in alphas") {
    Eigen::ArrayXd alphas(n);
    alphas << 0, 0, 1, 1, 1;
    Eigen::MatrixXd result = utils::populate_interaction_matrix(interaction, alphas, n);
    REQUIRE_FALSE(result.array().isFinite().all());
    result(0, 1) = 0;
    REQUIRE(result.array().isFinite().all());
  }
}

TEST_CASE("get_independent_row_indices; values", "[utils][matrix_utils]") {

  SECTION("all independent rows") {
    Eigen::MatrixXd m(3,5);
    m <<
      1, 0, 0, 1, 3,
      0, 1, 0, 1, 3,
      0, 0, 2, 0, 3;
    std::vector<Eigen::Index> expected = {0, 1, 2};
    std::vector<Eigen::Index> result = utils::get_independent_row_indices(m);
    REQUIRE(result == expected);
  }
  SECTION("some independent rows") {
    Eigen::MatrixXd m(4,4);
    m <<
      1, 0, 1, 3,
      0, 1, 1, 3,
      1, 0, 1, 3,
      0, 1, 1, 3;
    std::vector<Eigen::Index> expected = {0, 1};
    std::vector<Eigen::Index> result = utils::get_independent_row_indices(m);
    REQUIRE(result == expected);
  }
  SECTION("One independent row") {
    Eigen::MatrixXd m(2, 3);
    m <<
      2, 1, 4,
      2, 1, 4;
    std::vector<Eigen::Index> expected = {0};
    std::vector<Eigen::Index> result = utils::get_independent_row_indices(m);
    REQUIRE(result == expected);
  }
  SECTION("No independent rows - zero matrix") {
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(5, 5);
    std::vector<Eigen::Index> result = utils::get_independent_row_indices(m);
    REQUIRE(result.empty());
  }
}

TEST_CASE("complete_basis", "[utils][matrix_utils]") {
  SECTION("Basis already complete (square)") {
    Eigen::Index n = 3;
    Eigen::MatrixXd basis = Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd result = utils::complete_basis(basis);
    REQUIRE(result.rows() == n);
    REQUIRE(result.cols() == n);
    REQUIRE(result.isIdentity());
  }
  SECTION("One row basis") {
    Eigen::MatrixXd basis(1, 3);
    basis << -1, 0, 1;
    Eigen::MatrixXd expected(3, 3);
    expected <<
      -1, 0, 1,
       1, 0, 0,
       0, 1, 0;
    Eigen::MatrixXd result = utils::complete_basis(basis);
    REQUIRE(result.rows() == 3);
    REQUIRE(result.cols() == 3);
    REQUIRE((result.array() == expected.array()).all());
    Eigen::FullPivLU<Eigen::MatrixXd> lu(result);
    REQUIRE(lu.rank() == 3);
  }
  SECTION("two independent rows") {
    Eigen::MatrixXd basis(2, 3);
    basis <<
      1, 2, 3,
      0, 1, 4;
    Eigen::MatrixXd expected(3, 3);
    expected <<
      1, 2, 3,
      0, 1, 4,
      1, 0, 0;
    Eigen::MatrixXd result = utils::complete_basis(basis);
    REQUIRE(result.rows() == 3);
    REQUIRE(result.cols() == 3);
    REQUIRE((result.array() == expected.array()).all());
    Eigen::FullPivLU<Eigen::MatrixXd> lu(result);
    REQUIRE(lu.rank() == 3);
  }
  SECTION("Zero basis (empty matrix)") {
    Eigen::MatrixXd basis(0, 5);
    Eigen::MatrixXd result = utils::complete_basis(basis);
    REQUIRE(result.rows() == 5);
    REQUIRE(result.cols() == 5);
    REQUIRE(result.isIdentity());
  }
  // Only used from reaction_basis? What about square but not full rank.
}
