/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <cmath>
#include "burnman/eos/hp.hpp"
#include "burnman/eos/einstein.hpp"

double HP_TMT::compute_relative_thermal_pressure(
  double temperature,
  const MineralParams& params
) const {
  return compute_thermal_pressure(temperature, params)
    - compute_thermal_pressure(*params.T_0, params);
}

double HP_TMT::compute_thermal_pressure(
  double temperature,
  const MineralParams& params
) const {
  double E_th = einstein::compute_thermal_energy(
    temperature, *params.T_einstein, *params.napfu);
  double C_V0 = einstein::compute_molar_heat_capacity_v(
    *params.T_0, *params.T_einstein, *params.napfu);
  return (*params.a_0) * (*params.K_0) / C_V0 * E_th;
}
