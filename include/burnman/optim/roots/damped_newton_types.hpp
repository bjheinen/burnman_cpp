/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_OPTIM_ROOTS_DAMPED_NEWTON_TYPES_HPP_INCLUDED
#define BURNMAN_OPTIM_ROOTS_DAMPED_NEWTON_TYPES_HPP_INCLUDED

#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include "burnman/utils/constants.hpp"

namespace burnman {
namespace optim{
namespace roots{

/**
 * @brief Type for linear equality constraints (A·x + b <= 0)
 */
using LinearConstraints = std::pair<Eigen::MatrixXd, Eigen::VectorXd>;

/**
 * @brief Type for lambda bounda pair (min, max)
 */
using LambdaBounds = std::pair<double, double>;

/**
 * @brief Type for LambdaBounds function: returns (min_lambda, max_lambda).
 */
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

} // namespace roots
} // namespace optim
} // namespace burnman

#endif // BURNMAN_OPTIM_ROOTS_DAMPED_NEWTON_TYPES_HPP_INCLUDED
