/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_OPTIM_ROOTS_DAMPED_NEWTON_STATE_HPP_INCLUDED
#define BURNMAN_OPTIM_ROOTS_DAMPED_NEWTON_STATE_HPP_INCLUDED

#include <functional>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include "burnman/optim/roots/damped_newton_types.hpp"

namespace optim{
namespace roots{

/**
 * @struct DampedNewtonSolverState
 * @brief Stores internal DampedNewtonSolver state.
 */
struct DampedNewtonSolverState {
  const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& F_func;
  const LinearConstraints& linear_constraints;
  Eigen::VectorXd x;                        ///< Current solution vector.
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
  double dx_norm;                           ///< L2 norm of dx.
  double dxbar_j_norm;                      ///< L2 norm of dxbar_j.
  double h;                                 ///< Heuristic used to compute lambda.
  double lambda = 0.0;                      ///< Current step scaling (damping) factor.
  LambdaBounds lambda_bounds;               ///< Current bounds (min, max) on lambda.
  std::vector<std::pair<Eigen::Index, double>> violated_constraints;
  int n_iterations = 0;
  bool converged = false;
  bool minimum_lambda = false;
  bool persistent_bound_violation = false;
  bool require_posteriori_loop = true;
  // Constructor to store const refs
  DampedNewtonSolverState(
    const Eigen::VectorXd& x0,
    const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& F_func_,
    const LinearConstraints& linear_constraints_)
      : F_func(F_func_),
        linear_constraints(linear_constraints_),
        x(x0),
        F(F_func_(x0)),
        dxbar(Eigen::VectorXd::Ones(x0.size())),
        dx_prev(Eigen::VectorXd::Ones(x0.size())),
        n_constraints(linear_constraints_.first.rows())
    {}
};

} // namespace roots
} // namespace optim

#endif // BURNMAN_OPTIM_ROOTS_DAMPED_NEWTON_STATE_HPP_INCLUDED
