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
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "burnman/tools/equilibration/equality_constraints.hpp"
#include <Eigen/Dense>
#include "tolerances.hpp"
#include "solution_fixtures.hpp"

using namespace Catch::Matchers;
using namespace burnman;

TEST_CASE_METHOD(PyroliteAssemblageFixture, "Basic Constraints", "[tools][equilibration][equality_constraints]") {

  Eigen::VectorXd x(3);
  x << 5.0, 500.0, 2.0;
  Eigen::Index x_size = 3;
  Eigen::VectorXd deriv_ref = Eigen::VectorXd::Zero(x_size);

  SECTION("PressureConstraint") {
    auto c = equilibration::make_constraint<equilibration::PressureConstraint>(5.0);
    REQUIRE(c->evaluate(x, assemblage) == 0.0);
    deriv_ref(0) = 1.0;
    auto deriv_calc = c->derivative(x, assemblage, x_size);
    REQUIRE(deriv_calc.isApprox(deriv_ref, tol_rel));
  }
  SECTION("TemperatureConstraint") {
    auto c = equilibration::make_constraint<equilibration::TemperatureConstraint>(500.0);
    REQUIRE(c->evaluate(x, assemblage) == 0.0);
    deriv_ref(1) = 1.0;
    auto deriv_calc = c->derivative(x, assemblage, x_size);
    REQUIRE(deriv_calc.isApprox(deriv_ref, tol_rel));
  }
  SECTION("PTEllipseConstraint") {
    Eigen::Vector2d centre(5.0, 500.0);
    Eigen::Vector2d scaling(1.0, 5.0);
    double ref_eval = -1.0;
    auto c = equilibration::make_constraint<equilibration::PTEllipseConstraint>(centre, scaling);
    REQUIRE_THAT(c->evaluate(x, assemblage),
      WithinRel(ref_eval, tol_rel) || WithinAbs(ref_eval, tol_abs));
    auto test_data = GENERATE(
      std::make_pair(Eigen::Vector2d(1.0, 500.0), Eigen::Vector2d(-1.0, 0.0)),
      std::make_pair(Eigen::Vector2d(10.0, 500.0), Eigen::Vector2d(1.0, 0.0)),
      std::make_pair(Eigen::Vector2d(5.0, 100.0), Eigen::Vector2d(0.0, -0.2)),
      std::make_pair(Eigen::Vector2d(5.0, 1000.0), Eigen::Vector2d(0., 0.2))
    );
    x.segment(0, 2) = test_data.first;
    deriv_ref.segment(0, 2) = test_data.second;
    CAPTURE(x, deriv_ref);
    REQUIRE(c->derivative(x, assemblage, x_size).isApprox(deriv_ref, tol_rel));
  }
  SECTION("LinearXConstraint") {
    Eigen::VectorXd A(3);
    A << 0.2, 0.002, 0.5;
    double b = 3.0;
    auto c = equilibration::make_constraint<equilibration::LinearXConstraint>(A, b);
    REQUIRE(c->evaluate(x, assemblage) == 0);
    REQUIRE(c->derivative(x, assemblage, x_size).isApprox(A));
  }
}
