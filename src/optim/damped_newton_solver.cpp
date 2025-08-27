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
  // TODO: move assert of lambda bounds validity to function
  // at call of lambda_bounds_func
  // e.g. assert(lambda_bounds.second < 1.0 + eps);
  double lambda_j = std::min(1.0 / (h + this->settings.eps), lambda_bounds.second);
  return std::max(lambda_j, lambda_bounds.first);
}

bool DampedNewtonSolver::is_converged(
  const Eigen::VectorXd& dxbar_j,
  const Eigen::VectorXd& dx,
  double lambda,
  const LambdaBounds& lambda_bounds
) const {
  bool simplified_step_below_tol = (dxbar_j.array().abs() < this->settings.tol).all();
  bool full_step_below_tol = (dx.array().abs() < std::sqrt(10.0 * this->settings.tol)).all();
  bool lambda_is_max = std::abs(lambda - lambda_bounds.second) < this->settings.eps;
  return simplified_step_below_tol && full_step_below_tol && lambda_is_max;
}


std::vector<std::pair<int, double>>
DampedNewtonSolver::constrain_step_to_feasible_region(
  const Eigen::VectorXd& x,
  const Eigen::VectorXd& dx,
  double& lambda,
  Eigen::VectorXd& x_j
) {
  Eigen::VectorXd c_x = evaluate_constraints(x);
  Eigen::VectorXd c_x_j = evaluate_constraints(x_j);
  std::vector<std::pair<int, double>> violated_constraints;
  // TODO why separate n_constraints?
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(c_x.size()); ++i) {
    if (c_x_j(i) >= this->settings.eps) {
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

std::tuple<Eigen::VectorXd, Eigen::VectorXd, double>
DampedNewtonSolver::solve_subject_to_constraints(
  const Eigen::VectorXd& x,
  const Eigen::MatrixXd& Jx,
  const Eigen::VectorXd& c_x,
  const Eigen::MatrixXd& c_prime
) const {

  const Eigen::Index n_x = x.size();
  const Eigen::Index n_c = c_x.size();

  Eigen::MatrixXd JTJ_reg = Jx.transpose() * Jx
    + this->settings.regularisation * Eigen::MatrixXd::Identity(n_x, n_x);

  double norm = static_cast<double>(n_x * n_x) / JTJ_reg.norm();
  // KKT = np.block([[JTJ_reg * norm, c_prime.T], [c_prime, np.zeros((n_c, n_c))]])
  // Top left - JTJ_reg * norm
  // top right - c_prime.T
  // bottom left - c_prime
  // bottom right - zeros
  Eigen::MatrixXd KKT = Eigen::MatrixXd::Zero(n_x + n_c, n_x + n_c);
  KKT.topLeftCorner(n_x, n_x) = JTJ_reg * norm;
  KKT.topRightCorner(n_x, n_c) = c_prime.transpose();
  KKT.bottomLeftCorner(n_c, n_x) = c_prime;
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n_x + n_c);
  rhs.tail(n_c) = -c_x;
  // Get conditon number from SVD
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(KKT);
  double cond = svd.singularValues()(0)
    / svd.singularValues()(svd.singularValues().size() - 1);
  Eigen::VectorXd dx_lambda;
  if (cond < this->settings.condition_threshold_lu) {
    // LU solve
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(KKT);
    dx_lambda = lu.solve(rhs);
    // dx_lambda = KKT.partialPivLu().solve(rhs);
  } else if (cond < this->settings.condition_threshold_lstsq) {
    // LSTSQ solve in np
    // Here using the SVD seems closest equivalent
    dx_lambda = svd.solve(rhs);
  } else {
    // Python fallback is manual pseudo-inverse with SVD
    Eigen::VectorXd s = svd.singularValues();
    // cwiseInverse doesn't use a threshold
    Eigen::VectorXd s_inv = s.unaryExpr([](double val) {
      return (val > 1.0e-12) ? 1.0 / val : 0.0;
    });
    dx_lambda = svd.matrixV() * s_inv.asDiagonal()
      * svd.matrixU().transpose() * rhs;
    // This might just be equivalent to what svd.solve() does internally?
    // Other options - ColPivHouseholderQR
    // Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(KKT);
    // dx_lambda = qr.solve(rhs);
  }
  Eigen::VectorXd dx = dx_lambda.head(n_x);
  Eigen::VectorXd lambda = dx_lambda.tail(n_c);
  return {x + dx, lambdas, cond};
}

} // namespace roots
} // namespace optim
