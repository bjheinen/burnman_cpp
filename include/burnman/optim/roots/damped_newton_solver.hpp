/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_OPTIM_ROOTS_DAMPED_NEWTON_HPP_INCLUDED
#define BURNMAN_OPTIM_ROOTS_DAMPED_NEWTON_HPP_INCLUDED

#include <functional>
#include <tuple>
#include <utility>
#include <Eigen/Dense>
#include "burnman/optim/roots/damped_newton_types.hpp"
#include "burnman/optim/roots/damped_newton_state.hpp"

namespace optim{
namespace roots{

/**
 * @class DampedNewtonSolver
 * @brief Damped Newton solver with linear inequality constraints
 * TODO: Docs
 */
class DampedNewtonSolver {

 public:

  /**
   * @brief Construct solver with optional settings
   * @param solver_settings Solve parameters (default-intialised)
   */
  explicit DampedNewtonSolver(DampedNewtonSettings solver_settings = {})
    : settings(std::move(solver_settings)) {}

  /**
   * @brief Solve F(x) = 0 with optional linear constraints
   *
   * @param x0 Initial guess
   * @param F Objective function to solve
   * @param J Jacobian of F
   * @param linear_constraints Optional linear inequality constraints (A·x + b <= 0)
   *
   * @return DampedNewtonResult
   */
  [[nodiscard]] DampedNewtonResult solve(
    const Eigen::VectorXd& x0,
    const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& F,
    const std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& J,
    // TODO: default empty A/B or 0, -1?
    const LinearConstraints& linear_constraints = {
      Eigen::MatrixXd(0, 0),
      Eigen::VectorXd(0)
    }
  ) const;

 private:

  DampedNewtonSettings settings;

  /**
   * @brief Evaluates the linear constraints (A·x + b)
   */
  Eigen::VectorXd evaluate_constraints(
    const Eigen::VectorXd& x,
    const DampedNewtonSolverState& state
  ) const;

  /**
   * @brief Updates the damping factor, λ.
   *
   * Updates the damping factor for the current x and dx given a
   * heuristic, h, and the bounds lambda_bounds.
   *
   * Modifies state.lambda
   */
  void update_lambda(DampedNewtonSolverState& state) const;

  /**
   * @brief Solve a constrained Newton correction step using KKT system.
   *
   * Computes a step, `dx', that minimises the linearised residual ||J(x)·dx||
   * subject to linear equality constraints derived from the currently active
   * inequality constraints.
   *
   * The system is solved using the KKT (Karush-Kuhn-Tucker) formulation:
   * \f[
   *  \begin{bmatrix}
   *    J^T J + \alpha I & A^T \\
   *    A & 0
   *  \end{bmatrix}
   *  \begin{bmatrix}
   *    dx \\
   *    \lambda
   *  \end{bmatrix}
   *  =
   *  - \begin{bmatrix}
   *      0 \\
   *      c(x)
   *    \end{bmatrix}
   * \f]
   *
   * where:
   *  - \c J is the Jacobian at \c x
   *  - \c A is the constraint Jacobian (c_prime)
   *  - \c c(x) is the constraint evaluation
   *  - \c \lambda are the Lagrange multipliers
   *  - \c \alpha = \c settings.regularisation is an optional regularization parameter
   *
   * The KKT system is solved using one of thee strategies depending on
   * the estimated condition number of the matrix:
   *   1. TODO: decide on strategies in Eigen
   *
   * @param x Current solution vector.
   * @param Jx Current Jacobian matrix J(x).
   * @param c_x Values of active constraints at x.
   * @param c_prime Jacobian of active constraints (A in Ax + b = 0).
   *
   * @return tuple containing:
   *   x_new - Updated solution vector, x + dx.
   *   lambdas - Langrange multipliers for active constraints.
   *   condition_number - Estimated condition number of KKT matrix.
   */
  std::tuple<Eigen::VectorXd, Eigen::VectorXd, double> solve_subject_to_constraints(
    const Eigen::VectorXd& x,
    const Eigen::MatrixXd& Jx,
    const Eigen::VectorXd& c_x,
    const Eigen::MatrixXd& c_prime
  ) const;

  /**
   * @brief Project a trial Newton step back into feasible region
   *
   * The function checks if a trial point x_j = x + lambda·dx violates any
   * linear constraints (A·x + b <= 0). If so, it computes the maximum
   * allowable step scaling factor to remain feasible, reduces lambda
   * accordingly, and updates the trial iterate.
   *
   * The scaling factor is computed per violated constraint as
   *   lambda_i = c_x[i] / (c_x[i] - c_x_j[i])
   * where c_x and c_x_j are the constraint function values at x and x_j.
   * The smallest lambda_i is used to rescale the step to just touch the
   * first violated constraint.
   * @param state Current solver state. Uses:
   *   state.x Current solution vector
   *   state.dx Current Newton step direction
   *   state.lambda Current damping factor (scaled in place)
   *   state.x_j Current trial iterate, x + lambda·dx (updated in place)
   *   violated_constraints Sorted list given as (index, lambda_i) pairs,
   *     in ascending order of lambda_i.
   */
  void constrain_step_to_feasible_region(DampedNewtonSolverState& state) const;

  /**
   * @brief Attempts to find a constrained Newton step
   *
   * This function attempts to find a constrained Newton step when a step
   * along the standard Newton direction would violate active linear
   * inequality constraints (A·x + b <= 0).
   * Uses the method of Lagrange multipliers, attempting to "walk along"
   * the active constraints to remain in the feasible region while decreasing
   * the residual norm ||F(x)||.
   *
   * @param state Current solver state. Uses:
   *   state.x Current solution vector
   *   state.F Current function evaluation, F(x)
   *   state.dx_norm L2 norm of current Newton step
   *   state.luJ LU decomposition of the current Jacobian
   *   state.violated_constraints List of index, fraction pairs for constraints that would be violated by the current step.
   *   state.lambda Current damping factor (scaled in place).
   *   state.dx Current Newton step direction (modified in place).
   *   state.x_j Current trial iterate, x + lambda·dx (updated in place).
   *   state.persistent_bound_violation flag set.
   */
  void lagrangian_walk_along_constraints(DampedNewtonSolverState& state) const;

  /**
   * @brief Checks solve convergence
   *
   * Successful convergence requires meeting three criteria:
   *   - dx(simplified Newton step) < tol
   *   - dx(full Newton step) < sqrt(10*tol)
   *   - lambda = max (1 for full Newton step)
   */
  bool is_converged(const DampedNewtonSolverState& state) const;

  /**
   * @brief Posteriori step-size control loop
   *
   * Performs the posteriori step-size control loop of Deuflhard's
   * damped Newton method.
   *
   * After computing as trial step x_j = x + lambda·dx, this loop checks
   * whether the simplified Newton step (dx̄_j) decreases the residual norm
   * ||F(x)|| monotonically. If not, the damping factor lambda is reduced and
   * the trial step is recomputed until either monotonicity is restored
   * or the minimum lambda bound is reached. This procedure prevents divergence
   * and stabilised the Newton iteration in highly nonlinear regions.
   */
  void posteriori_loop(DampedNewtonSolverState& state) const;

  /**
   * @brief Sets solver termination info codes and message.
   */
  void make_termination_info(
    DampedNewtonResult& sol,
    const DampedNewtonSolverState& state
  ) const;

};

} // namespace roots
} // namespace optim

#endif // BURNMAN_OPTIM_ROOTS_DAMPED_NEWTON_HPP_INCLUDED
