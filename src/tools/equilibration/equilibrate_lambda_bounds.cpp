/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/tools/equilibration/equilibrate_lambda_bounds.hpp"
#include <cmath>
#include <algorithm>

std::pair<double, double> lambda_bounds_func(
  const Eigen::VectorXd& dx,
  const Eigen::VectorXd& x,
  const std::vector<int>& endmembers_per_phase
) {
  Eigen::ArrayXd max_steps = Eigen::ArrayXd::Constant(x.size(), 100000.0);
  // First two constraints are P & T - use biggest reasonable P-T steps
  max_steps(0) = 20.0e9;
  max_steps(1) = 500.0;
  int j = 2;
  for (int n : endmembers_per_phase) {
    if (x(j) + dx(j) < 0.0) {
      max_steps(j) = std::max(x(j)*0.999, 0.001);
    }
    for (int k = 1; k < n; ++k) {
      max_steps(j + k) = std::max(x(j + k) * 0.99, 0.01);
    }
    j += n;
  }
  double max_lambda = 1.0;
  for (Eigen::Index i = 0; i < dx.size(); ++i) {
    double step = std::abs(dx(i));
    double ratio = (step <= max_steps(i)) ? 1.0 : max_steps(i) / step;
    max_lambda = std::min(max_lambda, ratio);
  }
  return {1.0e-8, max_lambda};
}
