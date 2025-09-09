/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_SOLUTION_MODELS_IDEAL_INCLUDED
#define BURNMAN_CORE_SOLUTION_MODELS_IDEAL_HPP_INCLUDED

#include <Eigen/Dense>
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/core/solution_models/solution_model_base.hpp"

/**
 * @class IdealSolution
 * @brief Ideal solution model.
 *
 * Derived from SolutionModel.
 * Calculates the excess gibbs free energy and etropy due to configurational entropy.
 * Excess internal energy and volume are equal to zero.
 *
 * The multiplicity of each type of site in the structure is allowed to change
 * linearly as a function of endmember proportions. This makes the model equivalent
 * to the entropic part of a Temkin-type model (Temkin, 1945).
 *
 */
class IdealSolution : public SolutionModel{

 public:

  // Extend constructor
  IdealSolution(const PairedEndmemberList& endmember_list);

  // Public functions overriden from base class
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

  // Public or protected member variable?
  Eigen::ArrayXd endmember_configurational_entropies;
  Eigen::ArrayXd compute_endmember_configurational_entropies() const;

  // This is unused in python?
  double compute_configurational_entropy(
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::ArrayXd compute_ideal_excess_partial_gibbs(
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::ArrayXd compute_ideal_excess_partial_entropies(
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::ArrayXd compute_ideal_activities(
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::ArrayXd compute_log_ideal_activities(
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::MatrixXd compute_log_ideal_activity_derivatives(
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::MatrixXd compute_ideal_entropy_hessian(
    const Eigen::ArrayXd& molar_fractions) const;

  // Want to make and cache, ones, eyeones, eye
  // Not really used --> would make protected and create in setup

};

#endif // BURNMAN_CORE_SOLUTION_MODELS_IDEAL_HPP_INCLUDED
