/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/eos/models/bukowinski_electronic.hpp"
#include <cmath>

double bukowinski::compute_helmholtz_el(
  double temperature,
  double volume,
  const MineralParams& params
) {
  return -0.5 * (*params.bel_0)
    * std::pow(volume / *params.V_0, *params.gel)
    * (temperature * temperature - (*params.T_0) * (*params.T_0));
}

double bukowinski::compute_pressure_el(
  double temperature,
  double volume,
  const MineralParams& params
) {
  return 0.5 * (*params.gel) * (*params.bel_0)
    * std::pow(volume / *params.V_0, *params.gel)
    * (temperature * temperature - (*params.T_0) * (*params.T_0))
    / volume;
}

double bukowinski::compute_entropy_el(
  double temperature,
  double volume,
  const MineralParams& params
) {
  return *params.bel_0 * temperature
    * std::pow(volume / *params.V_0, *params.gel);
}

double bukowinski::compute_KT_over_V(
  double temperature,
  double volume,
  const MineralParams& params
) {
  return -(*params.gel - 1.0)
    * compute_pressure_el(temperature, volume, params)
    / volume;
}

double bukowinski::compute_CV_over_T(
  double volume,
  const MineralParams& params
) {
  return *params.bel_0 * std::pow(volume / *params.V_0, *params.gel);
}

double bukowinski::compute_alpha_KT(
  double temperature,
  double volume,
  const MineralParams& params
) {
  return *params.gel * (*params.bel_0) * temperature
    * std::pow(volume / *params.V_0, *params.gel)
    / volume;
}