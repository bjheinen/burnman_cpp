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
  const Eigen::VectorXd& [[maybe_unused]],
  const Assemblage& assemblage
) const {
  return assemblage.get_molar_entropy() * assemblage.get_n_moles() - this->value();
}

double VolumeConstraint::evaluate(
  const Eigen::VectorXd& [[maybe_unused]],
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

// Derivative implementations
Eigen::VectorXd PressureConstraint::derivative(
  const Eigen::VectorXd& x [[maybe_unused]],
  const Assemblage& assemblage [[maybe_unused]],
  Eigen::Index J_size
) const {
  Eigen::VectorXd row = Eigen::VectorXd::Zero(J_size);
  row(0) = 1.0;
  return row;
}

Eigen::VectorXd TemperatureConstraint::derivative(
  const Eigen::VectorXd& x [[maybe_unused]],
  const Assemblage& assemblage [[maybe_unused]],
  Eigen::Index J_size
) const {
  Eigen::VectorXd row = Eigen::VectorXd::Zero(J_size);
  row(1) = 1.0;
  return row;
}

Eigen::VectorXd EntropyConstraint::derivative(
  const Eigen::VectorXd& x,
  const Assemblage& assemblage,
  Eigen::Index J_size
) const {
  Eigen::VectorXd row = Eigen::VectorXd::Zero(J_size);
  row(0) = -assemblage.get_n_moles()
    * assemblage.get_thermal_expansivity()
    * assemblage.get_molar_volume();
  row(1) = assemblage.get_n_moles()
    * assemblage.get_molar_heat_capacity_p() / x(1);
  std::vector<int> embr_per_phase = assemblage.get_endmembers_per_phase();
  Eigen::Index j = 2;
  for (std::size_t k = 0; k < embr_per_phase.size(); ++k) {
    Eigen::Index n = static_cast<Eigen::Index>(embr_per_phase[k]);
    row(j) = assemblage.get_phase(k)->get_molar_entropy();
    if (n > 1) {
      auto ph = assemblage.get_phase<Solution>(k);
      // could do: auto ph  = std::static_pointer_cast<Solution>(phase); (with phase stored)
      row.segment(j + 1, n - 1) =
        assemblage.get_n_moles()
        * assemblage.get_molar_fractions()(k)
        * (ph->get_partial_entropies().tail(n - 1)
          - ph->get_partial_entropies()(0)).matrix();
    }
    j += n;
  }
  return row;
}

Eigen::VectorXd VolumeConstraint::derivative(
  const Eigen::VectorXd& x [[maybe_unused]],
  const Assemblage& assemblage,
  Eigen::Index J_size
) const {
  Eigen::VectorXd row = Eigen::VectorXd::Zero(J_size);
  row(0) = -assemblage.get_n_moles()
    * assemblage.get_molar_volume()
    / assemblage.get_isothermal_bulk_modulus_reuss();
  row(1) = assemblage.get_n_moles() * assemblage.get_molar_volume();
  std::vector<int> embr_per_phase = assemblage.get_endmembers_per_phase();
  Eigen::Index j = 2;
  for (std::size_t k = 0; k < embr_per_phase.size(); ++k) {
    Eigen::Index n = static_cast<Eigen::Index>(embr_per_phase[k]);
    row(j) = assemblage.get_phase(k)->get_molar_volume();
    if (n > 1) {
      auto ph = assemblage.get_phase<Solution>(k);
      row.segment(j + 1, n - 1) =
        assemblage.get_n_moles()
        * assemblage.get_molar_fractions()(k)
        * (ph->get_partial_volumes().tail(n-1)
          - ph->get_partial_volumes()(0)).matrix();
    }
    j += n;
  }
  return row;
}

Eigen::VectorXd PTEllipseConstraint::derivative(
  const Eigen::VectorXd& x,
  const Assemblage& assemblage [[maybe_unused]],
  Eigen::Index J_size
) const {
  Eigen::VectorXd row = Eigen::VectorXd::Zero(J_size);
  Eigen::Array2d v_scaled = (x.segment<2>(0) - centre).array() / scaling.array();
  row.segment<2>(0) = (v_scaled / (v_scaled.matrix().norm() * scaling.array())).matrix();
  return row;
}

Eigen::VectorXd LinearConstraintX::derivative(
  const Eigen::VectorXd& x [[maybe_unused]],
  const Assemblage& assemblage [[maybe_unused]],
  Eigen::Index J_size
) const {
  Eigen::VectorXd row = Eigen::VectorXd::Zero(J_size);
  row.head(A.size()) = A;
  return row;
}
