/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/core/solution.hpp"

void Solution::reset() {
  // Reset caches Material properties
  Material::reset();
  // Reset cached Solution properties
  excess_gibbs.reset();
  excess_volume.reset();
  excess_entropy.reset();
  excess_enthalpy.reset();
  activities.reset();
  activity_coefficients.reset();
  excess_partial_gibbs.reset();
  excess_partial_volumes.reset();
  excess_partial_entropies.reset();
  partial_gibbs.reset();
  partial_volumes.reset();
  partial_entropies.reset();
  gibbs_hessian.reset();
  entropy_hessian.reset();
  volume_hessian.reset();
}

double Solution::get_excess_gibbs() const {
  if (!excess_gibbs.has_value()) {
    excess_gibbs = compute_excess_gibbs();
  }
  return *excess_gibbs;
}

double Solution::get_excess_volume() const {
  if (!excess_volume.has_value()) {
    excess_volume = compute_excess_volume();
  }
  return *excess_volume;
}

double Solution::get_excess_entropy() const {
  if (!excess_entropy.has_value()) {
    excess_entropy = compute_excess_entropy();
  }
  return *excess_entropy;
}

double Solution::get_excess_enthalpy() const {
  if (!excess_enthalpy.has_value()) {
    excess_enthalpy = compute_excess_enthalpy();
  }
  return *excess_enthalpy;
}

Eigen::ArrayXd Solution::get_activities() const {
  if (!activities.has_value()) {
    activities = compute_activities();
  }
  return *activities;
}

Eigen::ArrayXd Solution::get_activity_coefficients() const {
  if (!activity_coefficients.has_value()) {
    activity_coefficients = compute_activity_coefficients();
  }
  return *activity_coefficients;
}

Eigen::ArrayXd Solution::get_excess_partial_gibbs() const {
  if (!excess_partial_gibbs.has_value()) {
    excess_partial_gibbs = compute_excess_partial_gibbs();
  }
  return *excess_partial_gibbs;
}

Eigen::ArrayXd Solution::get_excess_partial_volumes() const {
  if (!excess_partial_volumes.has_value()) {
    excess_partial_volumes = compute_excess_partial_volumes();
  }
  return *excess_partial_volumes;
}

Eigen::ArrayXd Solution::get_excess_partial_entropies() const {
  if (!excess_partial_entropies.has_value()) {
    excess_partial_entropies = compute_excess_partial_entropies();
  }
  return *excess_partial_entropies;
}

Eigen::ArrayXd Solution::get_partial_gibbs() const {
  if (!partial_gibbs.has_value()) {
    partial_gibbs = compute_partial_gibbs();
  }
  return *partial_gibbs;
}

Eigen::ArrayXd Solution::get_partial_volumes() const {
  if (!partial_volumes.has_value()) {
    partial_volumes = compute_partial_volumes();
  }
  return *partial_volumes;
}

Eigen::ArrayXd Solution::get_partial_entropies() const {
  if (!partial_entropies.has_value()) {
    partial_entropies = compute_partial_entropies();
  }
  return *partial_entropies;
}

Eigen::MatrixXd Solution::get_gibbs_hessian() const {
  if (!gibbs_hessian.has_value()) {
    gibbs_hessian = compute_gibbs_hessian();
  }
  return *gibbs_hessian;
}

Eigen::MatrixXd Solution::get_entropy_hessian() const {
  if (!entropy_hessian.has_value()) {
    entropy_hessian = compute_entropy_hessian();
  }
  return *entropy_hessian;
}

Eigen::MatrixXd Solution::get_volume_hessian() const {
  if (!volume_hessian.has_value()) {
    volume_hessian = compute_volume_hessian();
  }
  return *volume_hessian;
}
