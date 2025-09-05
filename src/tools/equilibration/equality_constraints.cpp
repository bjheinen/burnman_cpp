/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <numeric>
#include "burnman/tools/equilibrate/equality_constraints.hpp"
// TODO: could put prm in types.hpp, or a separate parameters.hpp
#include "burnman/tools/equilibrate/equilibration.hpp"

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

PhaseFractionConstraint::PhaseFractionConstraint(
  Eige::Index phase_index,
  double fraction,
  const EquilibrationParameters& prm)
  : LinearConstraintX(compute_A(phase_idx, fraction, prm), 0.0),
    phase_index(phase_index),
    fraction(fraction) {}

Eigen::VectorXd PhaseFractionConstraint::compute_A(
  Eigen::Index phase_index,
  double fraction,
  const EquilibrationParameters& prm
) {
  Eigen::VectorXd A = Eigen::VectorXd::Zero(prm.n_parameters);
  A(prm.phase_amount_indices) =
    Eigen::VectorXd::Constant(prm.phase_amount_indices.size(), -fraction);
  A(prm.phase_amount_indices(phase_idx)) += 1.0;
  return A;
}

PhaseCompositionConstraint::PhaseCompositionConstraint(
  Eigen::Index phase_index,
  const std::vector<std::string>& site_names,
  const Eigen::VectorXd& numerator,
  const Eigen::VectorXd& denominator,
  double value,
  const Assemblage& assemblage,
  const EquilibrationParameters& prm)
  : PhaseCompositionConstraint(
    compute_Ab(
      phase_index,
      site_names,
      numerator,
      denominator,
      value,
      assemblage,
      prm
    ),
    phase_index,
    site_names,
    numerator,
    denominator,
    value
  ) {}

PhaseCompositionConstraint::PhaseCompositionConstraint(
  const std::pair<Eigen::VectorXd, double>& Ab,
  Eigen::Index phase_index,
  const std::vector<std::string>& site_names,
  const Eigen::VectorXd& numerator,
  const Eigen::VectorXd& denominator,
  double value)
  : LinearConstraintX(Ab.first, Ab.second),
  phase_index(phase_index),
  site_names(site_names),
  numerator(numerator),
  denominator(denominator),
  value(value) {}

std::pair<Eigen::VectorXd, double> PhaseCompositionConstraint::compute_Ab(
  Eigen::Index phase_index,
  const std::vector<std::string>& site_names,
  const Eigen::VectorXd& numerator,
  const Eigen::VectorXd& denominator,
  double value,
  const Assemblage& assemblage,
  const EquilibrationParameters& prm
) {
  // Get phase
  const auto& phase = assemblage.get_phase<Solution>(static_cast<std::size_t>(phase_index));
  // If the phase isn't a solution we should throw an error
  if (!phase) {
    throw std::runtime_error(
      "Phase at index " + std::to_string(phase_index) + " is not a solution!");
  }
  // all solution site names
  const std::vector<std::string>& solution_site_names = phase.get_site_names();
  // Get site indices from site names
  Eigen::ArrayXi site_indices(static_cast<Eigen::Index>(site_names.size()));
  for (std::size_t i = 0; i < site_names.size(); ++i) {
    auto it = std::find(
      solution_site_names.begin(), solution_site_names.end(),
      site_names[i]
    );
    if (it == solution_site_names.end()) {
      throw std::runtime_error("Site name " + site_names[i] + " not found!")
    }
    site_indices(static_cast<Eigen::Index>(i)) =
      static_cast<int>(
        std::distance(solution_site_names.begin(), it)
      );
  }
  // Get indices for x vector
  const std::vector<int>& embr_per_phase = assemblage.get_endmembers_per_phase();
  Eigen::Index start_idx = std::accumulate(
    embr_per_phase.begin(),
    embr_per_phase.begin() + phase_index,
    0
  ) + 3;
  Eigen::Index n_indices = embr_per_phase[static_cast<std::size_t>(phase_index)] - 1;
  // Get endmember occupancy matrix
  const Eigen::ArrayXXd& n_occ = phase.get_endmember_n_occupancies();
  // Convert site constraints into endmember constraints
  // Shape should be (n_endmembers, 2)
  Eigen::MatrixXd numer_denom(2, numerator.size());
  numer_denom.row(0) = numerator;
  numer_denom.row(2) = denominator;
  Eigen::MatrixXd embr_c = n_occ(Eigen::all, site_indices).matrix()
    * numer_denom.transpose();
  // Compute b using 1st two values
  double b = value * embr_c(0, 1) - embr_c(0, 0);
  // Compute A
  embr_c.rowwise() -= embr_c.row(0).eval();
  Eigen::VectorXd numer = embr_c.col(0).segment(1, n_indices);
  Eigen::VectorXd denom = embr_c.col(1).segment(1, n_indices);
  Eigen::VectorXd A = Eigen::VectorXd::Zero(prm.n_parameters);
  A.segment(start_idx, n_indices) = numer - value * denom;
  return {A, b};
}

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
