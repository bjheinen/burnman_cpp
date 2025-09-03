/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_TOOLS_EQUILIBRATION_EQUALITY_CONSTRAINTS_HPP_INCLUDED
#define BURNMAN_TOOLS_EQUILIBRATION_EQUALITY_CONSTRAINTS_HPP_INCLUDED

#include <Eigen/Dense>
#include "burnman/core/assemblage.hpp"

// namespace?

// Virtual base class
/**
 * @brief Base class for linear equality constraints
 *
 * Use `make_constraint<ConstraintType>()' to construct a constraint.
 * Currently implemented constraints are:
 *   PressureConstraint
 *   TemperatureConstraint
 *   EntropyConstraint
 *   VolumeConstraint
 *   PTEllipseConstraint
 *   PhaseFractionConstraint
 *
 * Use `EqualityConstraint::evaluate()' to compute F.
 *
 */
class EqualityConstraint {
public:
  virtual ~EqualityConstraint() = default;
  virtual double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const = 0;
};

// Implemented constraints

class PressureConstraint : EqualityConstraint {
 public:
  explicit PressureConstraint(double value);
  double evaluate(const Eigen::VectorXd& x, const Assemblage&) const override;
 private:
  double value;
};

class TemperatureConstraint : EqualityConstraint {
 public:
  explicit TemperatureConstraint(double value);
  double evaluate(const Eigen::VectorXd& x, const Assemblage&) const override;
 private:
  double value;
};

class EntropyConstraint : EqualityConstraint {
 public:
  explicit EntropyConstraint(double value);
  double evaluate(const Eigen::VectorXd& x, const Assemblage&) const override;
 private:
  double value;
};

class VolumeConstraint : EqualityConstraint {
 public:
  explicit VolumeConstraint(double value);
  double evaluate(const Eigen::VectorXd& x, const Assemblage&) const override;
 private:
  double value;
};

class PTEllipseConstraint : EqualityConstraint {
 public:
  PTEllipseConstraint(const Eigen::Array2d& center, const Eigen::Vector2d& scaling);
  double evaluate(const Eigen::VectorXd& x, const Assemblage&) const override;
 private:
  Eigen::Vector2d center;
  Eigen::Vector2d scaling;
};

class LinearConstraintX : EqualityConstraint {
 public:
  LinearConstraintX(const Eigen::VectorXd& A, double b);
  double evaluate(const Eigen::VectorXd& x, const Assemblage&) const override;
 private:
  Eigen::VectorXd A;
  double b;
};

#endif // BURNMAN_TOOLS_EQUILIBRATION_EQUALITY_CONSTRAINTS_HPP_INCLUDED
