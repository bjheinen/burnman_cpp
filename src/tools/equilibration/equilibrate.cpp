/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <algorithm>
#include <cassert>
#include <cmath>
#include <Eigen/Dense>
#include "burnman/tools/equilibration/equilibrate.hpp"

Eigen::VectorXd get_parameter_vector(
  const Assemblage& assemblage,
  int n_free_compositional_vectors = 0
) {
  int n_params = assemblage.get_n_endmembers() + 2 + n_free_compositional_vectors;
  Eigen::VectorXd params = Eigen::VectorXd::Zero(n_params);
  Eigen::ArrayXd n_moles_per_phase = assemblage.n_moles * assemblage.molar_fractions;
  // check pressure & temperature are set
  if (!assemblage.has_state()) {
    throw std::runtime_error("You need to set_state before getting parameters");
  }
  // Set P & T as first two params
  params(0) = assemblage.pressure;
  params(1) = assemblage.temperature;
  std::vector<int> embr_per_phase = assemblage.get_endmembers_per_phase();
  Eigen::Index j = 2;
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(assemblage.get_n_phases()); ++i) {
    params(j) = n_moles_per_phase(i);
    if (auto ph = assemblage.get_phase<Solution>(static_cast<std::size_t>(i))) {
      Eigen::Index n_embr = static_cast<Eigen::Index>(embr_per_phase[i] - 1); // skip first embr
      params.segment(j + 1, n_embr) = ph.get_molar_fractions().tail(n_embr);
    }
    j += embr_per_phase[static_cast<std::size_t>(i)];
  }
  return params;
}

Eigen::VectorXd get_endmember_amounts(
  const Assemblage& assemblage
) {
  Eigen::ArrayXd phase_amounts = assemblage.get_n_moles * assemblage.get_molar_fractions();
  Eigen::VectorXd abs_amounts(assemblage.get_n_endmembers());
  std::vector<int> embr_per_phase = assemblage.get_endmembers_per_phase();
  Eigen::Index j = 0;
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(assemblage.get_n_phases()); ++i) {
    if (auto ph = assemblage.get_phase<Solution>(static_cast<std::size_t>(i))) {
      abs_amounts.segment(j, j + static_cast<Eigen::Index>(embr_per_phase(i))) =
        phase_amounts(i) * ph.get_molar_fractions();
    } else {
      abs_amounts(j) = phase_amounts(i);
    }
    j += embr_per_phase[static_cast<std::size_t>(i)];
  }
  return abs_amounts;
}

void set_composition_and_state_from_parameters(
  Assemblage& assemblage,
  const Eigen::VectorXd& parameters
) {
  // Set P & T (first two parameters)
  assemblage.set_state(parameters(0), parameters(1));
  Eigen::Index n_phases = static_cast<Eigen::Index>(assemblage.get_n_phases());
  Eigen::ArrayXd phase_amounts = Eigen::ArrayXd::Zero(n_phases);
  Eigen::Index i = 2;
  for (Eigen::Index phase_idx = 0; phase_idx < n_phases; ++phase_idx) {
    phase_amounts(phase_idx) = parameters(i);
    // TODO: get_phase take Eigen::Index?
    if (auto ph = assemblage.get_phase<Solution>(static_cast<std::size_t>(phase_idx))) {
      // TODO: make public get_n_endmembers - + use internally
      Eigen::Index n_mbrs = static_cast<Eigen::Index>(ph.compute_n_endmembers());
      Eigen::ArrayXd f = Eigen::ArrayXd::Zero(n_mbrs);
      f.segment(1, n_mbrs - 1) = parameters.segment(i + 1, n_mbrs - 1);
      f(0) = 1.0 - f.tail(n_mbrs - 1).sum();
      ph.set_composition(f);
      i += n_mbrs;
    } else {
      ++i;
    }
  }
  assert((phase_amounts > -1.0e-8).all());
  phase_amounts = phase_amounts.abs();
  assemblage.set_n_moles = phase_amounts.sum();
  assemblage.set_fractions(phase_amounts / assemblage.get_n_moles());
}

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
