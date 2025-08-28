/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_OPTIM_DAMPED_NEWTON_HPP_INCLUDED
#define BURNMAN_OPTIM_DAMPED_NEWTON_HPP_INCLUDED

#include <functional>
#include <optional>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "burnman/utils/constants.hpp"

// Do optim::roots::brent etc. instead
namespace optim{
namespace roots{

/**
 * @struct DampedNewtonSolverState
 * @brief Stores internal solver state.
 */
struct SolverState {
  Eigen::VectorXd x                         ///< Current solution vector.
  Eigen::VectorXd F;                        ///< Current function evaluation, F(x).
  Eigen::VectorXd dx;                       ///< Current Newton step direction.
  Eigen::VectorXd dxbar;                    ///< Simplified Newton step.
  Eigen::VectorXd dx_prev;                  ///< Newton step from previous iteration.
  Eigen::VectorXd x_j;                      ///< Trial iterate for solution vector.
  Eigen::VectorXd c_x_j;                    ///< Constraints evaluated at x_j.
  Eigen::VectorXd F_j;                      ///< Function evaluated at x_j.
  Eigen::VectorXd dxbar_j;                  ///< Simplified Newton step at trial iterate.
  Eigen::MatrixXd J;                        ///< Current Jacobian matrix J(x).
  Eigen::PartialPivLU<Eigen::MatrixXd> luJ; ///< LU decomposition of J.
  Eigen::Index n_constraints;               ///< Number of constraints.
  LambdaBounds lambda_bounds;               ///< Current bounds (min, max) on lambda.
  double dx_norm;                           ///< L2 norm of dx.
  double dxbar_j_norm;                      ///< L2 norm of dxbar_j.
  double h;                                 ///< Heuristic used to compute lambda.
  double lambda = 0.0;                      ///< Current step scaling (damping) factor.
  bool converged = false;
  bool minimum_lambda = false;
  bool persistent_bound_violation = false;
  bool require_posteriori_loop = true;
};

/**
 * @struct DampedNewtonResult
 * @brief Result of the damped Newton solver.
 *
 * This object stores the result of `DampedNewtonSolver::solve',
 * which attempts to solve F(x) = 0, subject to linear inequality
 * constraints A·x + b <= 0.
 *
 * If the solver fails to converge, it terminates with an informative
 * status code and message.
 *
 * The result object optionally stores iteration history of the solver.
 */
struct DampedNewtonResult {
  /**
   * @brief Final solution vector (x).
   */
  Eigen::VectorXd x;

  /**
   * @brief Function evaluation at the solution, F(x).
   */
  Eigen::VectorXd F;

  /**
   * @brief Jacobian matrix evaluated at the solution.
   */
  Eigen::MatrixXd J;

  /**
   * @brief Euclidean (L2) norm of F(x).
   */
  double F_norm;

  /**
   * @brief Number of Newton iterations.
   */
  int n_iterations = 0;

  /**
   * @brief True if solver converged successfully.
   */
  bool success = false;

  /**
   * @brief Termination code.
   *
   * 0 : Successful convergence
   * 1 : Failure (lambda reached its minimum bound)
   * 2 : Failure (descent direction violates constraints)
   * 3 : Failure (maximum iterations reached)
   */
  int code = -1;

  /**
   * @brief Human-readable termination message.
   */
  std::string message;

  /**
   * @brief Optional iteration history (present if store_iterates==true).
   */
  struct Iterates {
    std::vector<Eigen::VectorXd> x;
    std::vector<Eigen::VectorXd> F;
    std::vector<double> lambda;
  };
  std::optional<Iterates> iteration_history;

};

// Type for LambdaBounds function: returns (min_lambda, max_lambda)
using LambdaBounds = std::pair<double, double>;
using LambdaBoundsFunc = std::function<LambdaBounds(
  const Eigen::VectorXd&, const Eigen::VectorXd&)>;

/**
 * @struct DampedNewtonSettings
 * @brief Parameters controlling DampedNewtonSolver
 */
struct DampedNewtonSettings {
  bool store_iterates = false;              ///< Store iteration history (at each step)
  int max_iterations = 100;                 ///< Maximum number of Newton iterations
  double tol = 1.0e-6;                      ///< Convergence tolerance on |F|
  double regularisation = 0.0;              ///< Regularization parameter for the KKT system in Lagrangian solves
  double condition_threshold_lu = 1e12;     ///< Condition number below which LU decomposition is considered stable
  double condition_threshold_lstsq = 1e15;  ///< Condition number below which least-squares fallback is considered stable
  double eps = 2.0 * constants::precision::double_eps;
  LambdaBoundsFunc lambda_bounds_func =              ///< Callable (dx, x) that returns min, max for the damping factor
    [](const Eigen::VectorXd&, const Eigen::VectorXd&) {
      return std::make_pair(1.0e-8, 1.0);
    };
};

/**
 * @brief Type for linear equality constraints (A·x + b <= 0)
 */
using LinearConstraints = std::pair<Eigen::MatrixXd, Eigen::VectorXd>;

/**
 * @class DampedNewtonSolver
 * @brief Damped Newton solver with linear inequality constraints
 * TODO: Docs
 */
class DampedNewtonSolver {

 public:

  /**
   * @brief Construct solver with optional settings
   * @param settings Solve parameters (default-intialised)
   */
  explicit DampedNewtonSolver(DampedNewtonSettings settings = {})
    : settings(std::move(settings)) {}

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
    LinearConstraints linear_constraints = {
      Eigen::MatrixXd(0, 0),
      Eigen::VectorXd(0)
    }
  ) const;

 private:

  DampedNewtonSettings settings;
  LinearConstraints linear_constraints;

  /**
   * @brief Eavluates the linear constraints (A·x + b)
   */
  Eigen::VectorXd evaluate_constraints(const Eigen::VectorXd& x) const;

  /**
   * @brief Computes the damping factor, λ.
   *
   * Updates the damping factor for the current x and dx given a
   * heuristic, h, and the bounds lambda_bounds.
   */
  double compute_lambda(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& dx,
    double h,
    const LambdaBounds& lambda_bounds
  ) const;

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
   *
   * @param[in] x Current solution vector
   * @param[in] dx Current Newton step direction
   * @param[in,out] lambda Current damping factor (scaled in place)
   * @param[in,out] x_j Current trial iterate, x + lambda·dx (updated in place)
   *
   * @return violated_constraints Sorted list given as (index, lambda_i) pair,
   *         in ascending order of lambda_i.
   *
   * @note Modifies @p lambda and @p x_j.
   */
  std::vector<std::pair<int, double>>
  constrain_step_to_feasible_region(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& dx,
    double& lambda,
    Eigen::VectorXd& x_j
  );

  // TODO --> move x_j, lambda etc. to solver state sol?
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
   * @param[in] sol Current solver state with x and F
   * @param[in] dx_norm L2 norm of current Newton step
   * @param[in] luJ LU decomposition of the current Jacobian
   * @param[in] violated_constraints List of index, fraction pairs for constraints that would be violated by the current step.
   * @param[in,out] lambda Current damping factor (scaled in place).
   * @param[in,out] dx Current Newton step direction (modified in place).
   * @param[in,out] x_j Current trial iterate, x + lambda·dx (updated in place).
   *
   * @return persistent_bound_violation flag.
   */
  bool lagrangian_walk_along_constraints(
    const DampedNewtonResult& sol,
    const double dx_norm,
    const Eigen::PartialPivLU<Eigen::MatrixXd>& luJ,
    const std::vector<std::pair<int, double>>& violated_constraints
    double& lambda,
    Eigen::VectorXd& dx,
    Eigen::VectorXd& x_j
  );


  /**
   * @brief Checks solve convergence
   *
   * Successful convergence requires meeting three criteria:
   *   - dx(simplified Newton step) < tol
   *   - dx(full Newton step) < sqrt(10*tol)
   *   - lambda = max (1 for full Newton step)
   */
  bool is_converged(
    const Eigen::VectorXd& dxbar_j,
    const Eigen::VectorXd& dx,
    double lambda,
    const LambdaBounds& lambda_bounds
  ) const;

  // TODO --> think about return here (maybe just update state)
  std::tuple<?> posteriori_loop(
    Eigen::VectorXd& x,
    Eigen::VectorXd& Fx,
    Eigen::VectorXd& dx,
    double dx_norm,
    Eigen::VectorXd& dxbar_j,
    double dxbar_j_norm,
    Eigen::VectorXd& x_j,
    ? luJ,
    double lambda,
    LambdaBounds lambda_bounds,
    bool converged,
    bool minimum_lambda,
    bool persistent_bound_violation,
    bool require_posteriori_loop
  ) const;

  termination_info;
  // --> update state of Result object? or pass back termination info type and add to result?

};

} // namespace roots
} // namespace optim


#endif // BURNMAN_OPTIM_DAMPED_NEWTON_HPP_INCLUDED
