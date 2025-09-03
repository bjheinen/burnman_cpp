/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/tools/equilibrate/equality_constraints.hpp"

// TODO: think about namespaces

// Constraint constructors
PressureConstraint::PressureConstraint(double value)
  : value(value) {}

TemperatureConstraint::TemperatureConstraint(double value)
  : value(value) {}

EntropyConstraint::EntropyConstraint(double value)
  : value(value) {}

VolumeConstraint::VolumeConstraint(double value)
  : value(value) {}

PTEllipseConstraint::PTEllipseConstraint(
  const Eigen::Vector2d& centre,
  const Eigen::Vector2d& scaling)
  : centre(centre), scaling(scaling) {}

LinearConstraintX::LinearConstraintX(const Eigen::VectorXd& A, double b)
  : A(A), b(b) {}

// Evaluate implementations
double PressureConstraint::evaluate(
  const Eigen::VectorXd& x,
  const Assemblage& assemblage [[maybe_unused]]
) const {
  return x(0) - this->value;
}

double TemperatureConstraint::evaluate(
  const Eigen::VectorXd& x,
  const Assemblage& assemblage [[maybe_unused]]
) const {
  return x(1) - this->value;
}

double EntropyConstraint::evaluate(
  const Eigen::VectorXd&,
  const Assemblage& assemblage
) const {
  return assemblage.get_molar_entropy() * assemblage.get_n_moles() - this->value();
}

double VolumeConstraint::evaluate(
  const Eigen::VectorXd&,
  const Assemblage& assemblage
) const {
  return assemblage.get_molar_volum() * assemblage.get_n_moles() - this->value();
}

double PTEllipseConstraint::evaluate(
  const Eigen::VectorXd& x,
  const Assemblage& assemblage [[maybe_unused]]
) const {
  Eigen::Array2d v_scaled = (x.segment<2>(0).array() - this->centre.array()) / this->scaling.array();
  return v_scaled.matrix().norm() - 1.0;
}

double LinearConstraintX::evaluate(
  const Eigen::VectorXd& x,
  const Assemblage& assemblage [[maybe_unused]]
) const {
  return this->A.dot(x) - this->b;
}
