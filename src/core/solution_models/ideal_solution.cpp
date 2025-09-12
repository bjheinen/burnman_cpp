/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/core/solution_models/ideal_solution.hpp"
#include "burnman/utils/constants.hpp"
#include "burnman/utils/math_utils.hpp"

namespace burnman::solution_models {

// Constructor for IdealSolution
IdealSolution::IdealSolution(const PairedEndmemberList& endmember_list)
  : SolutionModel(endmember_list) {
  // Calculate configurational entropies also
  this->endmember_configurational_entropies = compute_endmember_configurational_entropies();
}

// Public function overrides for IdealSolution
Eigen::ArrayXd IdealSolution::compute_excess_partial_gibbs_free_energies(
  double pressure [[maybe_unused]],
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return compute_ideal_excess_partial_gibbs(temperature, molar_fractions);
}

Eigen::ArrayXd IdealSolution::compute_excess_partial_entropies(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  const Eigen::ArrayXd& molar_fractions
) const {
  return compute_ideal_excess_partial_entropies(molar_fractions);
}

Eigen::ArrayXd IdealSolution::compute_excess_partial_volumes(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  const Eigen::ArrayXd& molar_fractions [[maybe_unused]]
) const {
  return Eigen::ArrayXd::Zero(this->n_endmembers);
}

Eigen::MatrixXd IdealSolution::compute_gibbs_hessian(
  double pressure [[maybe_unused]],
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return -temperature * compute_ideal_entropy_hessian(molar_fractions);
}

Eigen::MatrixXd IdealSolution::compute_entropy_hessian(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  const Eigen::ArrayXd& molar_fractions
) const {
  return compute_ideal_entropy_hessian(molar_fractions);
}

Eigen::MatrixXd IdealSolution::compute_volume_hessian(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  const Eigen::ArrayXd& molar_fractions [[maybe_unused]]
) const {
  return Eigen::MatrixXd::Zero(this->n_endmembers, this->n_endmembers);
}

Eigen::ArrayXd IdealSolution::compute_activities(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  const Eigen::ArrayXd& molar_fractions
) const {
  return compute_ideal_activities(molar_fractions);
}

Eigen::ArrayXd IdealSolution::compute_activity_coefficients(
  double pressure [[maybe_unused]],
  double temperature [[maybe_unused]],
  const Eigen::ArrayXd& molar_fractions [[maybe_unused]]
) const {
  return Eigen::ArrayXd::Ones(this->n_endmembers);
}

// Private functions for IdealSolution
Eigen::ArrayXd IdealSolution::compute_endmember_configurational_entropies() const {
  return -constants::physics::gas_constant
  * (this->endmember_n_occupancies
    * (utils::logish(this->endmember_n_occupancies)
      - utils::logish(this->site_multiplicities))
    ).rowwise().sum();
}

Eigen::ArrayXd IdealSolution::compute_ideal_excess_partial_gibbs(
  double temperature,
  const Eigen::ArrayXd& molar_fractions
) const {
  return -(temperature * compute_ideal_excess_partial_entropies(molar_fractions));
}

Eigen::ArrayXd IdealSolution::compute_ideal_excess_partial_entropies(
  const Eigen::ArrayXd& molar_fractions
) const {
  return -(constants::physics::gas_constant
    * compute_log_ideal_activities(molar_fractions));
}

Eigen::ArrayXd IdealSolution::compute_ideal_activities(
  const Eigen::ArrayXd& molar_fractions
) const {
  // Dot product
  Eigen::ArrayXd reduced_n_occupancies = (this->endmember_n_occupancies.colwise() * molar_fractions).colwise().sum();
  Eigen::ArrayXd reduced_multiplicities = (this->site_multiplicities.colwise() * molar_fractions).colwise().sum();
  Eigen::ArrayXd reduced_occupancies = reduced_n_occupancies * utils::inverseish(reduced_multiplicities);
  Eigen::ArrayXd a = reduced_occupancies.transpose().replicate(this->n_endmembers, 1).pow(this->endmember_n_occupancies.array()).rowwise().prod();
  Eigen::ArrayXd norm_constants = (this->endmember_configurational_entropies / constants::physics::gas_constant).exp();
  return norm_constants * a;
}

Eigen::ArrayXd IdealSolution::compute_log_ideal_activities(
  const Eigen::ArrayXd& molar_fractions
) const {
  Eigen::ArrayXd reduced_n_occupancies = (this->endmember_n_occupancies.colwise() * molar_fractions).colwise().sum();
  Eigen::ArrayXd reduced_multiplicities = (this->site_multiplicities.colwise() * molar_fractions).colwise().sum();
  Eigen::ArrayXd lna = (
    this->endmember_n_occupancies.rowwise()
    * (utils::logish(reduced_n_occupancies)
      - utils::logish(reduced_multiplicities)
      ).transpose()
  ).rowwise().sum();
  Eigen::ArrayXd norm_constants = this->endmember_configurational_entropies / constants::physics::gas_constant;
  return lna + norm_constants;
}

Eigen::MatrixXd IdealSolution::compute_log_ideal_activity_derivatives(
  const Eigen::ArrayXd& molar_fractions
) const {
  Eigen::ArrayXd reduced_n_occupancies = (this->endmember_n_occupancies.colwise() * molar_fractions).colwise().sum();
  Eigen::ArrayXd reduced_multiplicities = (this->site_multiplicities.colwise() * molar_fractions).colwise().sum();
  Eigen::MatrixXd dlnadp =
    ((this->endmember_n_occupancies.rowwise() * utils::inverseish(reduced_n_occupancies).transpose()).matrix()
      * this->endmember_n_occupancies.matrix().transpose())
    - ((this->endmember_n_occupancies.rowwise() * utils::inverseish(reduced_multiplicities).transpose()).matrix()
      * this->site_multiplicities.matrix().transpose());
  return dlnadp;
}

Eigen::MatrixXd IdealSolution::compute_ideal_entropy_hessian(
  const Eigen::ArrayXd& molar_fractions
) const {
  return -(constants::physics::gas_constant
    * compute_log_ideal_activity_derivatives(molar_fractions));
}

} // namespace burnman::solution_models
