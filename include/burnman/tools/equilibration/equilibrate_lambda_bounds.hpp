/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_LAMBDA_BOUNDS_HPP_INCLUDED
#define BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_LAMBDA_BOUNDS_HPP_INCLUDED

#include <utility>
#include <vector>
#include <Eigen/Dense>

/**
 * @brief Returns lambda bounds for damped Newton solver.
 *
 * Computes bounds on lambda for the damped affine invariant modification to
 * Newton's method for nonlinear problems (Deuflhard, 1974;1975;2004).
 *
 * @param dx The proposed Newton step.
 * @param x The parameter vector.
 * @returns (min, max)
 */
std::pair<double, double> lambda_bounds_func(
  const Eigen::VectorXd& dx,
  const Eigen::VectorXd& x,
  const std::vector<int>& endmembers_per_phase);

#endif // BURNMAN_TOOLS_EQUILIBRATION_EQUILIBRATE_LAMBDA_BOUNDS_HPP_INCLUDED
