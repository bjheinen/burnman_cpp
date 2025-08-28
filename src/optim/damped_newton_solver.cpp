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

DampedNewtonResult DampedNewtonSolver::solve(
  const Eigen::VectorXd& x0,
  const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& F,
  const std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& J,
  LinearConstraints linear_constraints
) const {

  // Make solution object
  DampedNewtonResult sol;
  // Populate with starting vector etc.
  sol.x = x0;
  sol.F = F(x0);

  // TODO: factor out store_iteration?
  if (settings.store_iterates) {
    sol.iteration_history.x.push_back(sol.x);
    sol.iteration_history.F.push_back(sol.F);
    sol.iteration_history.lambda.push_back(0.0);
  }

  double lambda = 0.0;
  Eigen::VectorXd dxprev = Eigen::VectorXd::Ones(sol.x.size());
  Eigen::VectorXd dxbar = Eigen::VectorXd::Ones(sol.x.size());
  sol.n_iterations = 0;
  Eigen::Index n_constraints = evaluate_constraints(sol.x).size();
  bool minimum_lambda = False;
  bool converged = False;
  bool persistent_bound_violation = False;

  // TODO: perhaps move all solver info to sol --> then set as class member
  // Main loop
  while (sol.n_iterations < this->settings.max_iterations
      && !minimum_lambda
      && !persistent_bound_violation
      && !converged)
  {
    sol.J = J(sol.x);
    Eigen::PartialPivLU<Eigen::MatrixXd> luJ(sol.J);
    // TODO: declare this before? is ok just local to loop iteration?
    Eigen::VectorXd dx = luJ.solve(-sol.F);
    double dx_norm = dx.norm();
    // Get lambda_bounds from function
    LambdaBounds lambda_bounds = settings.lambda_bounds_func(dx, sol.x);
    double h = lambda
      * (dxbar - dx).norm()
      * dx_norm
      / (dxprev.norm * dxbar.norm);
    lambda = compute_lambda(sol.x, dx, h, lambda_bounds);

    Eigen::VectorXd x_j = sol.x + lambda * dx;
    Eigen::VectorXd c_x_j = evaluate_constraints(x_j);

    if ((c_x_j.array() >= this->settings.eps).any()) {
      // Updates lmabda and x_j in place
      std::vector<std::pair<int, double>> violated_constraints =
        constrain_step_to_feasible_region(
          sol.x, dx, lambda, x_j
        );
    }
    if (lambda < this->settings.eps) {
      TYPE? persistent_bound_violation =
        lagrangian_walk_along_constraints(
          sol, dx, luJ, dx_norm, violated_constraints
        );
    }

    Eigen::VectorXd F_j = F(x_j);
    Eigen::VectorXd dxbar_j = luJ.solve(-F_j);
    double dxbar_j_norm = dxbar_j.norm();
    converged = is_converged(dxbar_j, dx, lambda, lambda_bounds);
    bool require_posteriori_loop = !converged;

    // call posteriori_loop here

    // Store history
    // TODO: factor out store_iteration?
    if (this->settings.store_iterates) {
      sol.iteration_history.x.push_back(sol.x);
      sol.iteration_history.F.push_back(sol.F);
      sol.iteration_history.lambda.push_back(lambda);
    }
    // bump n_it
    ++sol.n_iterations;
  }

  // Final adjustments for constraints
  // TODO --> persist.. define outside of loop
  if (converged && !persistent_bound_violation) {
    sol.x = x_j + dxbar_j;
    if ((evaluate_constraints(sol.x).array() > 0.0).any()) {
      sol.x -= dxbar_j;
    }
  }
  sol.F = F(sol.x);
  sol.F_norm = sol.F.norm();
  sol.J = J(sol.x);

  // Call get_termination_info here

  return sol;
}

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

bool DampedNewtonSolver::lagrangian_walk_along_constraints(
  const DampedNewtonResult& sol,
  const double dx_norm,
  const Eigen::PartialPivLU<Eigen::MatrixXd>& luJ,
  const std::vector<std::pair<int, double>>& violated_constraints,
  double& lambda,
  Eigen::VectorXd& dx,
  Eigen::VectorXd& x_j
) {

  // TODO: use Eigen::Index for constraint indices?
  // Desparate active/inactive constraints
  std::vector<int> act_ind_temp;
  std::vector<int> inact_ind_temp;
  for (const auto& vc : violated_constraints) {
    int i = vc.first;
    double frac = vc.second;
    if (frac < this->settings.eps) {
      act_ind_temp.push_back(i);
    } else {
      inact_ind_temp.push_back(i);
    }
  }
  Eigen::VectorXi active_constraint_indices =
    Eigen::Map<Eigen::VectorXi>(act_ind_temp.data(), act_ind_temp.size());
  Eigen::VectorXi inactive_constraint_indices =
    Eigen::Map<Eigen::VectorXi>(inact_ind_temp.data(), inact_ind_temp.size());

  Eigen::VectorXd x_n = sol.x + dx;
  Eigen::VectorXd c_newton = evaluate_constraints(x_n)(active_constraint_indices);
  Eigen::MatrixXd c_A = (this->linear_constraints.first)(active_constraint_indices, Eigen::all);
  bool persistent_bound_violation = false;

  Eigen::VectorXd x_m;
  if (c_A.rows() > 0 && Eigen::FullPivLU<Eigen::MatrixXd>(c_A).rank() == dx.size()) {
    std::size_t n_act = active_constraint_indices.size();
    for (std::size_t i_rm = 0; i_rm < n_act; ++i_rm) {
      std::vector<int> p_act_ind_temp;
      for (std::size_t i = 0; i < n_act; ++i) {
        if (i != i_rm) {
          p_act_ind_temp.push_back(active_constraint_indices[i]);
        }
      }
      Eigen::VectorXi potential_active_indices = Eigen::Map<Eigen::VectorXi>(p_act_ind_temp.data(), p_act_ind_temp.size());
      c_newton = evaluate_constraints(sol.x + dx)(potential_active_indices);
      c_A = (this->linear_constraints.first)(potential_active_indices, Eigen::all);
      x_m = std::get<0>(solve_subject_to_constraints(x_n, sol.J, c_newton, c_A));
      if (evaluate_constraints(x_m)(active_constraint_indices[i_rm]) < 0.0) {
        break;
      }
    }
  } else {
    x_m = std::get<0>(solve_subject_to_constraints(x_n, sol.J, c_newton, c_A));
  }

  // Update dx, lambda
  dx = x_m - sol.x;
  LambdaBounds lambda_bounds_new = this->settings.lambda_bounds_func(dx, sol.x);
  lambda = lambda_bounds_new.second;
  x_j = sol.x + lambda * dx;
  // Check feasibility
  Eigen::VectorXd x_j_min = sol.x + lambda_bounds_new.first * dx;
  // TODO --> make F, J members?
  Eigen::VectorXd F_j_min = F(x_j_min);
  Eigen::VectorXd dxbar_j_min = luJ.solve(-F_j_min);
  double dxbar_j_min_norm = dxbar_j_min.norm();
  if (dxbar_j_min_norm > dx_norm || dx.norm() < this->settings.eps) {
    persistent_bound_violation = true;
  }

  // Check newly violated inactive constraints
  std::size_t n_inactive = inactive_constraint_indices.size();
  Eigen::VectorXd c_x_j = evaluate_constraints(x_j)(inactive_constraint_indices);
  if ((c_x_j.array() >= this->settings.eps).any()) {
    Eigen::VectorXd c_x = evaluate_constraints(sol.x)(inactive_constraint_indices);
    std::vector<std::pair<int, double>> newly_violated_constraints;
    for (std::size_t k = 0; k < n_inactive; ++k) {
      if (c_x_j[k] >= this->settings.eps) {
        double frac = c_x[k] / (c_x[k] - c_x_j[k]);
        newly_violated_constraints.emplace_back(inactive_constraint_indices[k], frac);
      }
    }
    std::sort(newly_violated_constraints.begin(), newly_violated_constraints.end(),
      [](auto& a, auto& b) { return a.second < b.second; });
    if (!newly_violated_constraints.empty()) {
      lambda *= newly_violated_constraints[0].second;
      x_j = sol.x + lambda * dx;
    }
  }
  return persistent_bound_violation;
}

} // namespace roots
} // namespace optim
