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
#include "burnman/eos/slb.hpp"
#include "burnman/eos/debye.hpp"
#include "burnman/eos/birch_murnaghan.hpp"
#include "burnman/optim/brent_solver.hpp"
#include "burnman/utils/constants.hpp"

bool SLB3::validate_parameters(MineralParams& params) {
  ;// TODO
}

double SLB3::compute_volume(
  double pressure,
  double temperature,
  const MineralParams& params
) const {
  ;
}

double SLB3::compute_pressure(
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

double SLB3::compute_grueneisen_parameter(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

double SLB3::compute_slb_grueneisen_parameter(
  double x,
  const MineralParams& params
) {
  ;
}

double SLB3::compute_isothermal_bulk_modulus_reuss(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

double SLB3::compute_isentropic_bulk_modulus_reuss(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

// Second order (SLB2)
double SLB2::compute_shear_modulus(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

// Third order (SLB3)
double SLB3::compute_shear_modulus(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

double SLB3::compute_molar_heat_capacity_v(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

double SLB3::compute_molar_heat_capacity_p(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

double SLB3::compute_thermal_expansivity(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

double SLB3::compute_gibbs_free_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

double SLB3::compute_entropy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

double SLB3::compute_molar_internal_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

double SLB3::compute_helmholtz_free_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

double SLB3::compute_enthalpy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  ;
}

double SLB3::compute_debye_temperature(
  double x,
  const MineralParams& params
) {
  ;
}

double compute_volume_dependent_q(
  double x,
  const MineralParams& params
) {
  ;
}

double compute_isotropic_eta_s(
  double x,
  const MineralParams& params
) {
  ;
}