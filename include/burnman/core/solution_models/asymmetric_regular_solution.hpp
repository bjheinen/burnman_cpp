/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_SOLUTION_MODELS_ASYMMETRIC_REGULAR_INCLUDED
#define BURNMAN_CORE_SOLUTION_MODELS_ASYMMETRIC_REGULAR_INCLUDED

#include <vector>
#include <Eigen/Dense>
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/core/solution_models/ideal_solution.hpp"

namespace burnman {
namespace solution_models {

/**
 * @class AsymmetricRegularSolution
 * @brief Asymmetric regular solution model.
 *
 * Derived from IdealSolution.
 * Implements the asymmetric regular solution model desribed in Holland &
 * Powell, 2003.
 *
 * The excess non-configurational Gibbs energy is given by the expression:
 * \f[
 * \mathcal{G}_{\textrm{excess}} = \alpha^T p (\phi^T W \phi)
 * \f]
 * \f$\alpha\f$ is a vector of van Laar parameters governing asymmetry in the
 * excess properties, i.e.:
 * \f[
 * \phi_i = \frac{\alpha_i p_i}{\sum_{k=1}^{n} \alpha_k p_k}, \quad
 * W_{ij} = \frac{2 w_{ij}}{\alpha_i + \alpha_j} \quad \text{for } i < j
 * \f]
 *
 * See `SymmetricRegularSolution' for the special case where all alpha = 1.
 */
class AsymmetricRegularSolution : public IdealSolution{

 public:

  AsymmetricRegularSolution(
    const types::PairedEndmemberList& endmember_list,
    std::vector<double> alphas,
    std::vector<std::vector<double>> energy_interaction,
    std::vector<std::vector<double>> volume_interaction = {},
    std::vector<std::vector<double>> entropy_interaction = {});

  // Alphas (van Laar) and interaction parameters
  Eigen::ArrayXd alphas;
  Eigen::MatrixXd W_e; // energy interactions
  Eigen::MatrixXd W_s; // entropy interactions
  Eigen::MatrixXd W_v; // volume interactions

  Eigen::ArrayXd compute_excess_partial_gibbs_free_energies(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::ArrayXd compute_excess_partial_entropies(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::ArrayXd compute_excess_partial_volumes(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::ArrayXd compute_activities(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::ArrayXd compute_activity_coefficients(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::MatrixXd compute_gibbs_hessian(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::MatrixXd compute_entropy_hessian(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::MatrixXd compute_volume_hessian(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

 private:

  Eigen::ArrayXd compute_phi(
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::ArrayXd compute_non_ideal_interactions(
    const Eigen::MatrixXd& W,
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::ArrayXd compute_non_ideal_excess_partial_gibbs(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::MatrixXd compute_non_ideal_hessian(
    const Eigen::MatrixXd& interactions,
    const Eigen::ArrayXd& molar_fractions) const;

};

} // namespace solution_models
} // namespace burnman

#endif // BURNMAN_CORE_SOLUTION_MODELS_ASYMMETRIC_REGULAR_INCLUDED
