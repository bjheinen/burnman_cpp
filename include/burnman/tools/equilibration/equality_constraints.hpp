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

#include <memory>
#include <utility>
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
 *   LinearConstraintX
 *
 * Use `EqualityConstraint::evaluate()' to compute F.
 * Use `EqualityConstraint::derivative()' to compute J.
 *
 */
class EqualityConstraint {
public:
  virtual ~EqualityConstraint() = default;
  virtual double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const = 0;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const = 0;
};

/**
 * @brief Helper function to construct constraint objects.
 */
template <typename ConstraintT, typename... Args>
std::unique_ptr<EqualityConstraint> make_constraint(Args&&... args) {
  return std::make_unique<ConstraintT>(std::forward<Args>(args)...);
}

/**
 * @brief Helper function to construct multiple constraints from vector of constraint values.
 */
template <typename ConstraintT>
std::vector<std::unique_ptr<EqualityConstraint>>
make_constraints_from_array(const Eigen::ArrayXd& values) {
  std::vector<std::unique_ptr<EqualityConstraint>> constraints;
  constraints.reserve(static_cast<std::size_t>(values.size()));
  for (double v : values) {
    constraints.push_back(make_constraint<ConstraintT>(v));
  }
  return constraints;
}

// Specialised template for PTEllipseConstraint
template <>
std::vector<std::unique_ptr<EqualityConstraint>>
make_constraints_from_array<PTEllipseConstraint>(
  const std::pair<Eigen::ArrayXd, Eigen::ArrayXd>& values
) {
  const Eigen::ArrayXd& centres = values.first;
  const Eigen::ArrayXd& scales = values.second;
  std::vector<std::unique_ptr<EqualityConstraint>> constraints;
  constraints.reserve(static_cast<std::size_t>(centres.size()));
  for (Eigen::Index i = 0; i < centres.szie(); ++i) {
    constraints.push_back(make_constraint<PTEllipseConstraint>(centres(i), scales(i)));
  }
  return constraints;
}

// Implemented constraints
class PressureConstraint : EqualityConstraint {
 public:
  explicit PressureConstraint(double value);
  double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const override;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const = 0;
 private:
  double value;
};

class TemperatureConstraint : EqualityConstraint {
 public:
  explicit TemperatureConstraint(double value);
  double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const override;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const = 0;
 private:
  double value;
};

class EntropyConstraint : EqualityConstraint {
 public:
  explicit EntropyConstraint(double value);
  double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const override;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const = 0;
 private:
  double value;
};

class VolumeConstraint : EqualityConstraint {
 public:
  explicit VolumeConstraint(double value);
  double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const override;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const = 0;
 private:
  double value;
};

class PTEllipseConstraint : EqualityConstraint {
 public:
  PTEllipseConstraint(const Eigen::Vector2d& centre, const Eigen::Vector2d& scaling);
  double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const override;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const = 0;
 private:
  Eigen::Vector2d centre;
  Eigen::Vector2d scaling;
};

class LinearConstraintX : EqualityConstraint {
 public:
  LinearConstraintX(const Eigen::VectorXd& A, double b);
  double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const override;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const = 0;
 private:
  Eigen::VectorXd A;
  double b;
};

// Subclasses for constructing Linear constraints from higher level args

// Forward declaration of EquilibrationParameters
struct EquilibrationParameters;

class PhaseFractionConstraint : LinearConstraintX {
 public:
  PhaseFractionConstraint(
    Eigen::Index phase_index,
    double fraction,
    const EquilibrationParameters& prm);
 private:
  Eigen::Index phase_index;
  double phase_fraction;
  static Eigen::VectorXd compute_A(
    std::size_t phase_idx,
    double fraction,
    const EquilibrationParameters& prm);
};

class PhaseCompositionConstraint : LinearConstraintX {
 public:
  PhaseCompositionConstraint(
    Eigen::Index phase_index,
    const std::vector<std::string>& site_names,
    const Eigen::VectorXd& numerator,
    const Eigen::VectorXd& denominator,
    double value,
    const Assemblage& assemblage,
    const EquilibrationParameters& prm);
 private:
  Eigen::Index phase_index;
  std::vector<std::string> site_names;
  Eigen::VectorXd numerator;
  Eigen::VectorXd denominator;
  double value;
  PhaseCompositionConstraint(
    const std::pair<Eigen::VectorXd, double>& Ab,
    Eigen::Index phase_index,
    const std::vector<std::string>& site_names,
    const Eigen::VectorXd& numerator,
    const Eigen::VectorXd& denominator,
    double value);
  static std::pair<Eigen::VectorXd, double> compute_Ab(
    Eigen::Index phase_index,
    const std::vector<std::string>& site_names,
    const Eigen::VectorXd& numerator,
    const Eigen::VectorXd& denominator,
    double value,
    const Assemblage& assemblage,
    const EquilibrationParameters& prm);
};

#endif // BURNMAN_TOOLS_EQUILIBRATION_EQUALITY_CONSTRAINTS_HPP_INCLUDED
