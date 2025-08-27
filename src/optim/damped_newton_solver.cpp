/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/optim/damped_newton_solver.hpp"
#include "burnman/utils/constants.hpp"

namespace optim{
namespace roots{

Eigen::VectorXd DampedNewtonSolver::evaluate_constraints(
  const Eigen::VectorXd& x
) const {
  const auto& [A, b] = this->linear_constraints;
  return A * x + b;
}

double DampedNewtonSolver::compute_lambda(
  const Eigen::VectorXd& x,
  const Eigen::VectorXd& dx,
  double h,
  const LambdaBounds& lambda_bounds
) const {
  // TODO: move eps to settings or constructor
  const double eps = 2 * constants::precision::double_eps;
  // TODO: move assert of lambda bounds validity to function
  // at call of lambda_bounds_func
  // e.g. assert(lambda_bounds.second < 1.0 + eps);
  double lambda_j = std::min(1.0 / (h + eps), lambda_bounds.second);
  return std::max(lambda_j, lambda_bounds.first);
}

bool DampedNewtonSolver::is_converged(
  const Eigen::VectorXd& dxbar_j,
  const Eigen::VectorXd& dx,
  double lambda,
  const LambdaBounds& lambda_bounds
) const {
  const double eps = 2.0 * constants::precision::double_eps; // TODO: from settings

  bool simplified_step_below_tol = (dxbar_j.array().abs() < this->settings.tol).all();
  bool full_step_below_tol = (dx.array().abs() < std::sqrt(10.0 * this->settings.tol)).all();
  bool lambda_is_max = std::abs(lambda - lambda_bounds.second) < eps;
  return simplified_step_below_tol && full_step_below_tol && lambda_is_max;
}


std::vector<std::pair<int, double>>
DampedNewtonSolver::constrain_step_to_feasible_region(
  const Eigen::VectorXd& x,
  const Eigen::VectorXd& dx,
  double& lambda,
  Eigen::VectorXd& x_j
) {
  // TODO grab from settings
  const double eps = 2.0 * constants::precision::double_eps;
  Eigen::VectorXd c_x = evaluate_constraints(x);
  Eigen::VectorXd c_x_j = evaluate_constraints(x_j);
  std::vector<std::pair<int, double>> violated_constraints;
  // TODO why separate n_constraints?
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(c_x.size()); ++i) {
    if (c_x_j(i) >= eps) {
      double lambda_i = c_x(i) / (c_x(i) - c_x_j(i));
      violated_constraints.emplace_back(i, lambda_i);
    }
  }
  if (!violated_constraints.empty()) {
    // Sort list
    std::sort(violated_constraints.begin(), violated_constraints.end(),
      [](const auto& a, const auto& b) {
        return a.second < b.second;
      }
    );
    // Update lambda and x_j
    lambda *= violated_constraints.front().second;
    x_j = x + lambda * dx;
  }
  // Return violated_constraints - lambda & x_j modified in place
  return violated_constraints;
}

} // namespace roots
} // namespace optim
