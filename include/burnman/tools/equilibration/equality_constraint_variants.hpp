/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_TOOLS_EQUILIBRATION_EQUALITY_CONSTRAINTS_VARIANTS_HPP_INCLUDED
#define BURNMAN_TOOLS_EQUILIBRATION_EQUALITY_CONSTRAINTS_VARIANTS_HPP_INCLUDED

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include "burnman/tools/equilibration/equality_constraint_base.hpp"
#include "burnman/tools/equilibration/equilibrate_types.hpp"

// Forward declarations
struct EquilibrationParameters;
class Assemblage;

class PressureConstraint : EqualityConstraint {
 public:
  explicit PressureConstraint(double value);
  std::unique_ptr<EqualityConstraint> clone() const override;
  double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const override;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const override;
 protected:
  double value;
};

class TemperatureConstraint : EqualityConstraint {
 public:
  explicit TemperatureConstraint(double value);
  std::unique_ptr<EqualityConstraint> clone() const override;
  double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const override;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const override;
 protected:
  double value;
};

class EntropyConstraint : EqualityConstraint {
 public:
  explicit EntropyConstraint(double value);
  std::unique_ptr<EqualityConstraint> clone() const override;
  double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const override;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const override;
 protected:
  double value;
};

class VolumeConstraint : EqualityConstraint {
 public:
  explicit VolumeConstraint(double value);
  std::unique_ptr<EqualityConstraint> clone() const override;
  double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const override;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const override;
 protected:
  double value;
};

class PTEllipseConstraint : EqualityConstraint {
 public:
  PTEllipseConstraint(const Eigen::Vector2d& centre, const Eigen::Vector2d& scaling);
  std::unique_ptr<EqualityConstraint> clone() const override;
  double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const override;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const override;
 protected:
  Eigen::Vector2d centre;
  Eigen::Vector2d scaling;
};

class LinearXConstraint : EqualityConstraint {
 public:
  LinearXConstraint(const Eigen::VectorXd& A, double b);
  std::unique_ptr<EqualityConstraint> clone() const override;
  double evaluate(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage) const override;
  virtual Eigen::VectorXd derivative(
    const Eigen::VectorXd& x,
    const Assemblage& assemblage,
    Eigen::Index J_size) const override;
 protected:
  Eigen::VectorXd A;
  double b;
};

// Higher level constraints that generate generic LinearXConstraints from
// user-friendly/readable args

class PhaseFractionConstraint : LinearXConstraint {
 public:
  PhaseFractionConstraint(
    Eigen::Index phase_index,
    double fraction,
    const EquilibrationParameters& prm);
  std::unique_ptr<EqualityConstraint> clone() const override;
 protected:
  Eigen::Index phase_index;
  double phase_fraction;
 private:
  static Eigen::VectorXd compute_A(
    std::size_t phase_idx,
    double fraction,
    const EquilibrationParameters& prm);
};

class PhaseCompositionConstraint : LinearXConstraint {
 public:
  PhaseCompositionConstraint(
    Eigen::Index phase_index,
    const std::vector<std::string>& site_names,
    const Eigen::VectorXd& numerator,
    const Eigen::VectorXd& denominator,
    double value,
    const Assemblage& assemblage,
    const EquilibrationParameters& prm);
  std::unique_ptr<EqualityConstraint> clone() const override;
 protected:
  Eigen::Index phase_index;
  std::vector<std::string> site_names;
  Eigen::VectorXd numerator;
  Eigen::VectorXd denominator;
  double value;
 private:
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

#endif // BURNMAN_TOOLS_EQUILIBRATION_EQUALITY_CONSTRAINTS_VARIANTS_HPP_INCLUDED
