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

class PressureConstraint : public EqualityConstraint {
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
  double get_value() const { return value_; };
 protected:
  double value_;
};

class TemperatureConstraint : public EqualityConstraint {
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
  double get_value() const { return value_; }
 protected:
  double value_;
};

class EntropyConstraint : public EqualityConstraint {
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
  double value_;
};

class VolumeConstraint : public EqualityConstraint {
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
  double value_;
};

class PTEllipseConstraint : public EqualityConstraint {
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
  Eigen::Vector2d get_scaling() const { return scaling_; }
 protected:
  Eigen::Vector2d centre_;
  Eigen::Vector2d scaling_;
};

class LinearXConstraint : public EqualityConstraint {
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
  Eigen::VectorXd A_;
  double b_;
};

// Higher level constraints that generate generic LinearXConstraints from
// user-friendly/readable args

class PhaseFractionConstraint : public LinearXConstraint {
 public:
  PhaseFractionConstraint(
    Eigen::Index phase_index,
    double phase_fraction,
    const EquilibrationParameters& prm);
  std::unique_ptr<EqualityConstraint> clone() const override;
 protected:
  Eigen::Index phase_index_;
  double phase_fraction_;
 private:
  static Eigen::VectorXd compute_A(
    Eigen::Index phase_index,
    double fraction,
    const EquilibrationParameters& prm);
};

class PhaseCompositionConstraint : public LinearXConstraint {
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
  Eigen::Index phase_index_;
  std::vector<std::string> site_names_;
  Eigen::VectorXd numerator_;
  Eigen::VectorXd denominator_;
  double value_;
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
