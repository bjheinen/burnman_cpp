/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
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
