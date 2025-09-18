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
#include <utility>
#include <Eigen/Dense>
#include "burnman/tools/equilibration/equilibrate_types.hpp"
#include "tolerances.hpp"
#include "solution_fixtures.hpp"

using namespace Catch::Matchers;
using namespace burnman;

TEST_CASE("Basic Constraints", "[tools][equilibration][equality_constraints]") {
  Assemblage assemblage;
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

TEST_CASE_METHOD(PyroliteAssemblageFixture, "Assemblage based constraints", "[tools][equilibration][equality_constraints]") {
  assemblage.set_state(50.e9, 1000.0);
  assemblage.set_n_moles(10);
  assemblage.set_method(types::EOSType::Auto);
  // Quick check that fixture hasn't changed
  REQUIRE(assemblage.get_n_endmembers() == 6);
  Eigen::VectorXd x = Eigen::VectorXd::Ones(8);
  Eigen::VectorXd deriv_ref(8);
  Eigen::Index x_size = 8;
  SECTION("EntropyConstraint") {
    double target_entropy = assemblage.get_molar_entropy() * assemblage.get_n_moles();
    auto c = equilibration::make_constraint<equilibration::EntropyConstraint>(target_entropy);
    REQUIRE(c->evaluate(x, assemblage) == 0);
    deriv_ref << -3.46442683e-09, 1.06218738e+03, 1.63441377e+02, 1.88155532e+02, 3.58896932e+02, 7.47098832e+01, 9.18054116e+01, 1.66273849e+02;
    // Relax tolerance, can make strict by using tol_rel instead
    REQUIRE(c->derivative(x, assemblage, x_size).isApprox(deriv_ref, 1.0e-07));
  }
  SECTION("VolumeConstraint") {
    double target_volume = assemblage.get_molar_volume() * assemblage.get_n_moles();
    auto c = equilibration::make_constraint<equilibration::VolumeConstraint>(target_volume);
    REQUIRE(c->evaluate(x, assemblage) == 0);
    deriv_ref << -1.01238991e-16, 2.00192406e-01, 2.14204549e-05, 6.72869416e-06, 3.44588368e-06, 1.00009341e-01, 2.00000226e+00, 2.37812313e-05;
    // Relax tolerance, can make strict by using tol_rel instead
    REQUIRE(c->derivative(x, assemblage, x_size).isApprox(deriv_ref, 1.0e-9));
  }
}

TEST_CASE_METHOD(PyroliteAssemblageFixture, "PhaseFractionConstraint", "[tools][equilibration][equality_constraints]") {
  assemblage.set_state(50.e9, 1000.0);
  assemblage.set_n_moles(10);
  assemblage.set_method(types::EOSType::Auto);
  types::FormulaMap composition = {{"Fe", 0.2}, {"Mg", 2.0}, {"Si", 1.9}, {"Ca", 0.2}, {"Al", 0.4}, {"O", 6.8}};
  // Quick check that fixture hasn't changed
  REQUIRE(assemblage.get_n_endmembers() == 6);
  Eigen::VectorXd x = Eigen::VectorXd::Ones(8);
  Eigen::Index x_size = 8;
  Eigen::VectorXd deriv_ref(8);
  equilibration::EquilibrationParameters prm; // = equilibration::get_equilibration_parameters(assemblage, composition, {});
  prm.n_parameters = 8;
  prm.phase_amount_indices = Eigen::ArrayXi(3);
  prm.phase_amount_indices << 2, 5, 7;
  SECTION("PhaseFractionConstraint") {
    // {bdg, fper, capv} == {0.7, 0.2, 0.1}
    Eigen::Index phase_index = 1;
    double phase_fraction = 0.2;
    auto c = equilibration::make_constraint<equilibration::PhaseFractionConstraint>(phase_index, phase_fraction, prm);
    double eval_ref = 0.4;
    REQUIRE_THAT(c->evaluate(x, assemblage),
      WithinRel(eval_ref, tol_rel) || WithinAbs(eval_ref, tol_abs));
    deriv_ref << 0, 0, -0.2, 0, 0, 0.8, 0, -0.2;
    REQUIRE(c->derivative(x, assemblage, x_size).isApprox(deriv_ref, tol_rel));
  }
  SECTION("PhaseCompositionConstraint") {
    Eigen::Index phase_index = 1;
    std::vector<std::string> site_names = {"Mg_A", "Fe_A"};
    Eigen::VectorXd numerator(2);
    Eigen::VectorXd denominator(2);
    numerator << 0, 1;
    denominator << 1, 1;
    double value = 0.1;
    auto c = equilibration::make_constraint<equilibration::PhaseCompositionConstraint>(
      phase_index, site_names, numerator, denominator, value, assemblage, prm);
    double eval_ref = 0.9;
    REQUIRE_THAT(c->evaluate(x, assemblage),
      WithinRel(eval_ref, tol_rel) || WithinAbs(eval_ref, tol_abs));
    deriv_ref << 0, 0, 0, 0, 0, 0, 1, 0;
    REQUIRE(c->derivative(x, assemblage, x_size).isApprox(deriv_ref, tol_rel));
  }
}

TEST_CASE("Test ConstraintGroup and List", "[tools][equilibration][equality_constraints]") {
  // Single constraints
  auto c_p = equilibration::make_constraint<equilibration::PressureConstraint>(5.0);
  auto c_t = equilibration::make_constraint<equilibration::TemperatureConstraint>(500.0);
  // Single value constraints
  Eigen::ArrayXd values = Eigen::ArrayXd::Random(10);
  equilibration::ConstraintGroup cg;
  REQUIRE_NOTHROW(cg = equilibration::make_constraints_from_array<equilibration::PressureConstraint>(values));
  REQUIRE_NOTHROW(cg = equilibration::make_constraints_from_array<equilibration::TemperatureConstraint>(values));
  REQUIRE_NOTHROW(cg = equilibration::make_constraints_from_array<equilibration::VolumeConstraint>(values));
  REQUIRE_NOTHROW(cg = equilibration::make_constraints_from_array<equilibration::EntropyConstraint>(values));
  // PTEllipse
  Eigen::Array<double, 2, Eigen::Dynamic> centres = Eigen::Array<double, 2, Eigen::Dynamic>::Random(2, 10);
  Eigen::Array<double, 2, Eigen::Dynamic> scales = Eigen::Array<double, 2, Eigen::Dynamic>::Random(2, 10);
  REQUIRE_NOTHROW(cg = equilibration::make_constraints_from_array<equilibration::PTEllipseConstraint>(std::make_pair(centres, scales)));
  // Currently no array constructor for LinearXConstraint and derived
  // Make list from previously constructed constraints and constraint groups
  equilibration::ConstraintList cl = equilibration::make_constraint_list(c_p, c_t, cg);
  REQUIRE(cl.size() == 3);
  REQUIRE(cl[2].size() == 10);
  REQUIRE_NOTHROW(
    cl = equilibration::make_constraint_list(
      equilibration::make_constraint<equilibration::PressureConstraint>(50.e9),
      equilibration::make_constraints_from_array<equilibration::TemperatureConstraint>(values),
      equilibration::make_constraints_from_array<equilibration::PTEllipseConstraint>(std::make_pair(centres,scales)),
      c_t
    )
  );
  REQUIRE(cl.size() == 4);
}
