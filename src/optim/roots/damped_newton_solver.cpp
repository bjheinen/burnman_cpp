/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/optim/roots/damped_newton_solver.hpp"
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <string>
#include <vector>

namespace optim{
namespace roots{

DampedNewtonResult DampedNewtonSolver::solve(
  const Eigen::VectorXd& x0,
  const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& F_func,
  const std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& J_func,
  const LinearConstraints& linear_constraints
) const {

  // Make solution result and solver iteration state objects
  DampedNewtonResult sol;
  DampedNewtonSolverState state{x0, F_func, linear_constraints};

  if (settings.store_iterates) {
    // Initialise Iterates in-place
    sol.iteration_history.emplace();
    sol.iteration_history->x.push_back(state.x);
    sol.iteration_history->F.push_back(state.F);
    sol.iteration_history->lambda.push_back(state.lambda);
  }

  // Main loop
  while (state.n_iterations < this->settings.max_iterations
      && !state.minimum_lambda
      && !state.persistent_bound_violation
      && !state.converged)
  {
    state.J = J_func(state.x);
    state.luJ.compute(state.J);
    state.dx = state.luJ.solve(-state.F);
    state.dx_norm = state.dx.norm();
    // Get lambda_bounds from function
    state.lambda_bounds = settings.lambda_bounds_func(state.dx, state.x);
    state.h = state.lambda
      * (state.dxbar - state.dx).norm()
      * state.dx_norm
      / (state.dx_prev.norm() * state.dxbar.norm());
    update_lambda(state);

    state.x_j = state.x + state.lambda * state.dx;
    state.c_x_j = evaluate_constraints(state.x_j, state);

    if ((state.c_x_j.array() >= this->settings.eps).any()) {
      // Updates lambda, x_j, violated_constraints in place
      constrain_step_to_feasible_region(state);
    }
    if (state.lambda < this->settings.eps) {
      // Updates dx, lambda, x_j, persistent_bound_violation
      lagrangian_walk_along_constraints(state);
    }

    state.F_j = F_func(state.x_j);
    state.dxbar_j = state.luJ.solve(-state.F_j);
    state.dxbar_j_norm = state.dxbar_j.norm();
    state.converged = is_converged(state);
    state.require_posteriori_loop = !state.converged;
    posteriori_loop(state);

    // Store history
    // TODO: factor out store_iteration?
    if (this->settings.store_iterates) {
      sol.iteration_history->x.push_back(state.x);
      sol.iteration_history->F.push_back(state.F);
      sol.iteration_history->lambda.push_back(state.lambda);
    }
    // bump n_it
    ++state.n_iterations;
  }

  // Final adjustments for constraints
  if (state.converged && !state.persistent_bound_violation) {
    state.x = state.x_j + state.dxbar_j;
    if ((evaluate_constraints(state.x, state).array() > 0.0).any()) {
      state.x -= state.dxbar_j;
    }
  }

  // Store result
  sol.x = state.x;
  sol.F = F_func(sol.x);
  sol.J = J_func(sol.x);
  sol.F_norm = sol.F.norm();
  sol.n_iterations = state.n_iterations;
  make_termination_info(sol, state);

  return sol;
}

Eigen::VectorXd DampedNewtonSolver::evaluate_constraints(
  const Eigen::VectorXd& x,
  const DampedNewtonSolverState& state
) const {
  const auto& [A, b] = state.linear_constraints;
  return A * x + b;
}

void DampedNewtonSolver::update_lambda(DampedNewtonSolverState& state) const {
  // TODO: move assert of lambda bounds validity to function
  // at call of lambda_bounds_func
  // e.g. assert(lambda_bounds.second < 1.0 + eps);
  double lambda_j = std::min(1.0 / (state.h + this->settings.eps), state.lambda_bounds.second);
  state.lambda = std::max(lambda_j, state.lambda_bounds.first);
}

bool DampedNewtonSolver::is_converged(
  const DampedNewtonSolverState& state
) const {
  bool simplified_step_below_tol = (state.dxbar_j.array().abs() < this->settings.tol).all();
  bool full_step_below_tol = (state.dx.array().abs() < std::sqrt(10.0 * this->settings.tol)).all();
  bool lambda_is_max = std::abs(state.lambda - state.lambda_bounds.second) < this->settings.eps;
  return simplified_step_below_tol && full_step_below_tol && lambda_is_max;
}

void DampedNewtonSolver::constrain_step_to_feasible_region(
  DampedNewtonSolverState& state
) const {
  Eigen::VectorXd c_x = evaluate_constraints(state.x, state);
  Eigen::VectorXd c_x_j = evaluate_constraints(state.x_j, state);
  state.violated_constraints.clear();
  // TODO: just use c_x.size()?
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(state.n_constraints); ++i) {
    if (c_x_j(i) >= this->settings.eps) {
      double lambda_i = c_x(i) / (c_x(i) - c_x_j(i));
      state.violated_constraints.emplace_back(i, lambda_i);
    }
  }
  if (!state.violated_constraints.empty()) {
    // Sort list
    std::sort(state.violated_constraints.begin(), state.violated_constraints.end(),
      [](const auto& a, const auto& b) {
        return a.second < b.second;
      }
    );
    // Update lambda and x_j
    state.lambda *= state.violated_constraints.front().second;
    state.x_j = state.x + state.lambda * state.dx;
  }
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
  Eigen::VectorXd lambdas = dx_lambda.tail(n_c);
  return {x + dx, lambdas, cond};
}

void DampedNewtonSolver::lagrangian_walk_along_constraints(
  DampedNewtonSolverState& state
) const {
  // Separate active/inactive constraints
  std::vector<int> act_ind_temp;
  std::vector<int> inact_ind_temp;
  for (const auto& vc : state.violated_constraints) {
    int i = static_cast<int>(vc.first);
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
  Eigen::VectorXd x_n = state.x + state.dx;
  Eigen::VectorXd c_newton = evaluate_constraints(x_n, state)(active_constraint_indices);
  Eigen::MatrixXd c_A = (state.linear_constraints.first)(active_constraint_indices, Eigen::all);
  state.persistent_bound_violation = false;
  Eigen::VectorXd x_m;
  if (c_A.rows() > 0 && Eigen::FullPivLU<Eigen::MatrixXd>(c_A).rank() == state.dx.size()) {
    Eigen::Index n_act = active_constraint_indices.size();
    for (Eigen::Index i_rm = 0; i_rm < n_act; ++i_rm) {
      std::vector<int> p_act_ind_temp;
      for (Eigen::Index i = 0; i < n_act; ++i) {
        if (i != i_rm) {
          p_act_ind_temp.push_back(active_constraint_indices[i]);
        }
      }
      Eigen::VectorXi potential_active_indices = Eigen::Map<Eigen::VectorXi>(p_act_ind_temp.data(), p_act_ind_temp.size());
      c_newton = evaluate_constraints(state.x + state.dx, state)(potential_active_indices);
      c_A = (state.linear_constraints.first)(potential_active_indices, Eigen::all);
      x_m = std::get<0>(solve_subject_to_constraints(x_n, state.J, c_newton, c_A));
      if (evaluate_constraints(x_m, state)(active_constraint_indices[i_rm]) < 0.0) {
        break;
      }
    }
  } else {
    x_m = std::get<0>(solve_subject_to_constraints(x_n, state.J, c_newton, c_A));
  }
  // Update dx, lambda
  state.dx = x_m - state.x;
  LambdaBounds lambda_bounds_new = this->settings.lambda_bounds_func(state.dx, state.x);
  state.lambda = lambda_bounds_new.second;
  state.x_j = state.x + state.lambda * state.dx;
  // Check feasibility
  Eigen::VectorXd x_j_min = state.x + lambda_bounds_new.first * state.dx;
  Eigen::VectorXd F_j_min = state.F_func(x_j_min);
  Eigen::VectorXd dxbar_j_min = state.luJ.solve(-F_j_min);
  double dxbar_j_min_norm = dxbar_j_min.norm();
  if (dxbar_j_min_norm > state.dx_norm || state.dx.norm() < this->settings.eps) {
    state.persistent_bound_violation = true;
  }
  // Check newly violated inactive constraints
  Eigen::Index n_inactive = inactive_constraint_indices.size();
  Eigen::VectorXd c_x_j = evaluate_constraints(state.x_j, state)(inactive_constraint_indices);
  if ((c_x_j.array() >= this->settings.eps).any()) {
    Eigen::VectorXd c_x = evaluate_constraints(state.x, state)(inactive_constraint_indices);
    std::vector<std::pair<Eigen::Index, double>> newly_violated_constraints;
    for (Eigen::Index k = 0; k < n_inactive; ++k) {
      if (c_x_j[k] >= this->settings.eps) {
        double frac = c_x[k] / (c_x[k] - c_x_j[k]);
        newly_violated_constraints.emplace_back(inactive_constraint_indices[k], frac);
      }
    }
    std::sort(newly_violated_constraints.begin(), newly_violated_constraints.end(),
      [](auto& a, auto& b) { return a.second < b.second; });
    if (!newly_violated_constraints.empty()) {
      state.lambda *= newly_violated_constraints[0].second;
      state.x_j = state.x + state.lambda * state.dx;
    }
  }
}

void DampedNewtonSolver::posteriori_loop(DampedNewtonSolverState& state) const {
  // Make local variables
  Eigen::VectorXd x_j = state.x_j;
  Eigen::VectorXd dxbar_j = state.dxbar_j;
  double dxbar_j_norm = state.dxbar_j_norm;
  while (
    state.require_posteriori_loop &&
    !state.minimum_lambda &&
    !state.persistent_bound_violation
  ) {
    if (dxbar_j_norm <= state.dx_norm) {
      if (dxbar_j_norm <= this->settings.eps) {
        state.converged = true;
      }
      state.x = x_j;
      state.F = state.F_func(x_j);
      state.dxbar = dxbar_j;
      state.dx_prev = state.dx;
      state.require_posteriori_loop = false;
    } else {
      if (std::abs(state.lambda - state.lambda_bounds.first) < this->settings.eps) {
        state.minimum_lambda = true;
      }
      double h_j = (2.0 / state.lambda) * (dxbar_j - (1.0 - state.lambda) * state.dx).norm() / state.dx_norm;
      double lambda_j = std::min(state.lambda_bounds.second, 1.0 / h_j);
      state.lambda = std::min(lambda_j, state.lambda / 2.0);
      state.lambda = std::max(state.lambda, state.lambda_bounds.first);
      x_j = state.x + state.lambda * state.dx;
      Eigen::VectorXd F_j = state.F_func(x_j);
      dxbar_j = state.luJ.solve(-F_j);
      dxbar_j_norm = dxbar_j.norm();
    }
  }
}

void DampedNewtonSolver::make_termination_info(
  DampedNewtonResult& sol,
  const DampedNewtonSolverState& state
) const {
  sol.success = state.converged;
  if (sol.success) {
    sol.code = 0;
    sol.message =
      "The solver successfully found a root after "
      + std::to_string(state.n_iterations)
      + " iterations.";
  } else {
    if (state.minimum_lambda) {
      sol.code = 1;
      sol.message =
        "The function is too non-linear for lower lambda bound ("
        + std::to_string(state.lambda_bounds.first)
        + ").";
    } else if (state.persistent_bound_violation) {
      sol.code = 2;
      std::string s = "[";
      for (std::size_t i = 0; i < state.violated_constraints.size(); ++i) {
        s += std::to_string(state.violated_constraints[i].first);
        if (i != state.violated_constraints.size() - 1) {
          s += ", ";
        }
      }
      s += "]";
      sol.message =
        std::string("The descent vector crosses the constraints with ")
        + "the following indices: "
        + s + "."
    } else if (state.n_iterations == this->settings.max_iterations) {
      sol.code = 3;
      sol.message =
        "The solver reached max_iterations ("
        + std::to_string(this->settings.max_iterations)
        + ").";
    } else {
      // Silently might be better? Return code=-1
      sol.code = -1;
      sol.message = "Error: Unknown termination of solver.";
    }
  }
}

} // namespace roots
} // namespace optim
