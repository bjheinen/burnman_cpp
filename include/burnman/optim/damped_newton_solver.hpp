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

#include <Eigen/Dense>
#include <functional>
#include <optional>
#include <string>
#include <vector>

// Do optim::roots::brent etc. instead
namespace optim{
namespace roots{

/**
 * @struct DampedNewtonResult
 * @brief Result of the damped Newton solver.
 *
 * This object stores the result of `DampedNewtonSolver::solve',
 * which attempts to solve F(x) = 0, subject to linea inequality
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
  LambdaBoundsFunc lambda_bounds =              ///< Callable (dx, x) that returns min, max for the damping factor
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

  Eigen::VectorXd constraints(const Eigen::VectorXd& x) const;

  double update_lambda(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& dx,
    const Eigen::VectorXd& h,
    LambdaBounds lambda_bounds
  ) const;

  std::tuple<Eigen::VectorXd, Eigen::VectorXd, double> solve_subject_to_constraints(
    const Eigen::VectorXd& x,
    const Eigen::MatrixXd& Jx,
    const Eigen::VectorXd& c_x,
    const Eigen::MatrixXd& c_prime
  ) const;

  std::tuple<double, Eigen::VectorXd, std::vector<std::pair<int, double>>>
  constrain_step_to_feasible_region(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& dx,
    int n_constraints,
    double lambda,
    Eigen::VectorXd& x_j
  ) const;

  std::pair<Eigen::VectorXd, double> lagrangian_walk_along_constraints(
    const DampedNewtonResult& sol,
    const Eigen::VectorXd& dx,
    double dx_norm,
    ? luJ,
    const std::vector<std::pair<int, double>>& violated_constraints
  ) const;

  bool is_converged(
    Eigen::VectorXd& dxbar_j,
    Eigen::VectorXd& dx,
    double lambda,
    LambdaBounds lambda_bounds
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
