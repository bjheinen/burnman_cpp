/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/core/solution_models/asymmetric_regular_solution.hpp"
#include <stdexcept>
#include "burnman/utils/constants.hpp"
#include "burnman/utils/matrix_utils.hpp"

namespace burnman::solution_models {

// Constructor for AsymmetricRegularSolution
AsymmetricRegularSolution::AsymmetricRegularSolution(
  const types::PairedEndmemberList& endmember_list,
  std::vector<double> alphas_vector,
  std::vector<std::vector<double>> energy_interaction,
  std::vector<std::vector<double>> volume_interaction,
  std::vector<std::vector<double>> entropy_interaction
) : IdealSolution(endmember_list) {
  // Map alphas to an Eigen::ArrayXd
  this->alphas = Eigen::Map<Eigen::ArrayXd>(
    alphas_vector.data(), alphas_vector.size());

  this->W_e = utils::populate_interaction_matrix(
    utils::jagged2square(energy_interaction, this->n_endmembers),
    this->alphas,
    this->n_endmembers);

  if (!volume_interaction.empty()) {
    this->W_v = utils::populate_interaction_matrix(
      utils::jagged2square(volume_interaction, this->n_endmembers),
      this->alphas,
      this->n_endmembers);
  } else {
    this->W_v = Eigen::MatrixXd::Zero(this->n_endmembers, this->n_endmembers);
  }
  if (!entropy_interaction.empty()) {
    this->W_s = utils::populate_interaction_matrix(
      utils::jagged2square(entropy_interaction, this->n_endmembers),
      this->alphas,
      this->n_endmembers);
  } else {
    this->W_s = Eigen::MatrixXd::Zero(this->n_endmembers, this->n_endmembers);
  }
}

// Public function overrides for AsymmetricRegularSolution
Eigen::ArrayXd AsymmetricRegularSolution::compute_excess_partial_gibbs_free_energies(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return IdealSolution::compute_excess_partial_gibbs_free_energies(
      pressure, temperature, molar_fractions)
    + compute_non_ideal_excess_partial_gibbs(
      pressure, temperature, molar_fractions);
}

Eigen::ArrayXd AsymmetricRegularSolution::compute_excess_partial_entropies(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return IdealSolution::compute_excess_partial_entropies(
      pressure, temperature, molar_fractions)
    + compute_non_ideal_interactions(this->W_s, molar_fractions);
}

Eigen::ArrayXd AsymmetricRegularSolution::compute_excess_partial_volumes(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return IdealSolution::compute_excess_partial_volumes(
      pressure, temperature, molar_fractions)
    + compute_non_ideal_interactions(this->W_v, molar_fractions);
}

Eigen::MatrixXd AsymmetricRegularSolution::compute_gibbs_hessian(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  Eigen::MatrixXd interactions = this->W_e - temperature * this->W_s + pressure * this->W_v;
  return IdealSolution::compute_gibbs_hessian(
      pressure, temperature, molar_fractions)
    + compute_non_ideal_hessian(interactions, molar_fractions);
}

Eigen::MatrixXd AsymmetricRegularSolution::compute_entropy_hessian(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return IdealSolution::compute_entropy_hessian(
      pressure, temperature, molar_fractions)
    + compute_non_ideal_hessian(this->W_s, molar_fractions);
}

Eigen::MatrixXd AsymmetricRegularSolution::compute_volume_hessian(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return IdealSolution::compute_volume_hessian(
      pressure, temperature, molar_fractions)
    + compute_non_ideal_hessian(this->W_v, molar_fractions);
}

Eigen::ArrayXd AsymmetricRegularSolution::compute_activities(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return IdealSolution::compute_activities(
      pressure, temperature, molar_fractions)
    * compute_activity_coefficients(
      pressure, temperature, molar_fractions);
}

Eigen::ArrayXd AsymmetricRegularSolution::compute_activity_coefficients(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  if (temperature < constants::precision::abs_tolerance) {
    throw std::runtime_error("Activity coefficients undefined at 0 K.");
  }
  return (
    compute_non_ideal_excess_partial_gibbs(
      pressure, temperature, molar_fractions)
    / (constants::physics::gas_constant * temperature)
  ).exp();
}

// Private compute functions for AsymmetricRegularSolution

Eigen::ArrayXd AsymmetricRegularSolution::compute_phi(
  const Eigen::ArrayXd& molar_fractions
) const {
  Eigen::ArrayXd phi = this->alphas * molar_fractions;
  return phi / phi.sum();
}

Eigen::ArrayXd AsymmetricRegularSolution::compute_non_ideal_interactions(
  const Eigen::MatrixXd& W,
  const Eigen::ArrayXd& molar_fractions
) const {
  Eigen::ArrayXd phi = compute_phi(molar_fractions);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(this->n_endmembers, this->n_endmembers);
  Eigen::MatrixXd q = I.rowwise() - phi.matrix().transpose();
  Eigen::ArrayXd W_int = -this->alphas * (q * W).cwiseProduct(q).rowwise().sum().array();
  return W_int;
}

Eigen::ArrayXd AsymmetricRegularSolution::compute_non_ideal_excess_partial_gibbs(
  double pressure,
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  Eigen::ArrayXd E_int = compute_non_ideal_interactions(this->W_e, molar_fractions);
  Eigen::ArrayXd S_int = compute_non_ideal_interactions(this->W_s, molar_fractions);
  Eigen::ArrayXd V_int = compute_non_ideal_interactions(this->W_v, molar_fractions);
  return E_int - temperature * S_int + pressure * V_int;
}

Eigen::MatrixXd AsymmetricRegularSolution::compute_non_ideal_hessian(
  const Eigen::MatrixXd& interactions,
  const Eigen::ArrayXd& molar_fractions
) const {
  // Maybe factor out - reused in non_ideal_interactions
  Eigen::ArrayXd phi = compute_phi(molar_fractions);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(this->n_endmembers, this->n_endmembers);
  Eigen::MatrixXd q = I.rowwise() - phi.matrix().transpose();
  double sum_pa = molar_fractions.matrix().dot(this->alphas.matrix());
  Eigen::MatrixXd alpha_outer_product = this->alphas.matrix() * (this->alphas/sum_pa).matrix().transpose();
  Eigen::MatrixXd qWq_product = q * interactions * q.transpose();
  Eigen::MatrixXd weighted_product = qWq_product.cwiseProduct(alpha_outer_product);
  Eigen::MatrixXd hessian = weighted_product + weighted_product.transpose();
  return hessian;
}

} // namespace burnman::solution_models
