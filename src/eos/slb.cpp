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
#include <stdexcept>
#include "burnman/eos/slb.hpp"
#include "burnman/eos/debye.hpp"
#include "burnman/eos/birch_murnaghan.hpp"
#include "burnman/optim/brent_solver.hpp"
#include "burnman/utils/constants.hpp"

bool SLB3::validate_parameters(MineralParams& params) {
  ;// TODO
}

std::pair<double, double> SLB3::get_b_g_el(
  const MineralParams& params [[maybe_unused]]
) const {
  return {0.0, 1.0};
}

std::pair<double, double> SLB3Conductive::get_b_g_el(
  const MineralParams& params
) const {
  return {*params.bel_0, *params.gel};
}

double SLB3::compute_volume(
  double pressure,
  double temperature,
  const MineralParams& params
) const {

  double gamma_0 = *params.grueneisen_0;
  double V_lo = *params.V_0 * 0.01;
  // Eq. 47
  double a1_ii = 6.0 * gamma_0;
  double a2_iikk = -12.0 * gamma_0
    + 36.0 * gamma_0 * gamma_0
    - 18.0 * (*params.q_0) * gamma_0;
  double b_iikk = 9.0 * (*params.K_0); // Eq.28
  double b_iikkmm = 27.0 * (*params.K_0) * (*params.Kprime_0 - 4.0); // Eq.29
  auto [bel_0, gel] = get_b_g_el(params);

  // TODO Bracketing




}

double SLB3::compute_pressure(
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double x = *params.V_0 / volume;
  double debye_temperature = compute_debye_temperature(x, params);
  double gamma = compute_slb_grueneisen_parameter(x, params);
  double E_th = debye::compute_thermal_energy(
    temperature, debye_temperature, *params.napfu);
  double E_th_ref = debye::compute_thermal_energy(
    *params.T_0, debye_temperature, *params.napfu);
  double b_iikk = 9.0 * (*params.K_0); // Eq.28
  double b_iikkmm = 27.0 * (*params.K_0) * (*params.Kprime_0 - 4.0); // Eq.29
  double x_cbrt = std::cbrt(x);
  double x_23 = x_cbrt * x_cbrt;
  double f = 0.5 * (x_23 - 1.0); // Eq.24
  double two_f_plus1 = 2.0 * f + 1;
  double two_f_plus1_52 = std::sqrt(two_f_plus1) * two_f_plus1 * two_f_plus1;
  constexpr double ONE_THIRD = 1.0 / 3.0;
  // Eq. 21
  return ONE_THIRD * two_f_plus1_52
    * ((b_iikk * f) + (0.5 * b_iikkmm * f * f))
    + gamma * (E_th - E_th_ref) / volume;
}

double SLB3::compute_grueneisen_parameter(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  // Simple alias to slb_g so SLB3Conductive can easily override
  double x = *params.V_0 / volume;
  return compute_slb_grueneisen_parameter(x, params);
}

double SLB3::compute_slb_grueneisen_parameter(
  double x,
  const MineralParams& params
) {
  // TODO:: Factor out Eq. 47 etc. to re-use (eta, q, etc.)
  double gamma_0 = *params.grueneisen_0;
  double x_cbrt = std::cbrt(x);
  double x_23 = x_cbrt * x_cbrt;
  double f = 0.5 * (x_23 - 1.0);
  // Eq. 47
  double a1_ii = 6.0 * gamma_0;
  double a2_iikk = -12.0 * gamma_0
    + 36.0 * gamma_0 * gamma_0
    - 18.0 * (*params.q_0) * gamma_0;
  // Eq. 41
  double nu_o_nu0_sq = 1.0 + a1_ii * f + 0.5 * a2_iikk * f * f;
  constexpr double ONE_SIXTH = 1.0 / 6.0;
  return ONE_SIXTH / nu_o_nu0_sq * (2.0 * f + 1.0) * (a1_ii + a2_iikk * f);
}

double SLB3::compute_isothermal_bulk_modulus_reuss(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double x = *params.V_0 / volume;
  double debye_temperature = compute_debye_temperature(x, params);
  double gamma = compute_slb_grueneisen_parameter(x, params);
  double q = compute_volume_dependent_q(x, params);
  double E_th = debye::compute_thermal_energy(
    temperature, debye_temperature, *params.napfu);
  double E_th_ref = debye::compute_thermal_energy(
    *params.T_0, debye_temperature, *params.napfu);
  double C_v = debye::compute_molar_heat_capacity_v(
    temperature, debye_temperature, *params.napfu);
  double C_v_ref = debye::compute_molar_heat_capacity_v(
    *params.T_0, debye_temperature, *params.napfu);
  return BM3::compute_bm_bulk_modulus(volume, params)
    + (gamma + 1.0 - q) * (gamma / volume) * (E_th - E_th_ref)
    - (gamma * gamma / volume)
      * (C_v * temperature - C_v_ref * (*params.T_0));
}

double SLB3::compute_isentropic_bulk_modulus_reuss(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double K_T = compute_isothermal_bulk_modulus_reuss(
    pressure,
    temperature,
    volume,
    params
  );
  double alpha = compute_thermal_expansivity(
    pressure,
    temperature,
    volume,
    params
  );
  double gamma = compute_grueneisen_parameter(
    pressure,
    temperature,
    volume,
    params);
  return K_T * (1.0 + gamma * alpha * temperature);
}

double SLB3::compute_shear_modulus_delta(
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double x = *params.V_0 / volume;
  double debye_temperature = compute_debye_temperature(x, params);
  double eta_s = compute_isotropic_eta_s(x, params);
  double E_th = debye::compute_thermal_energy(
    temperature, debye_temperature, *params.napfu);
  double E_th_ref = debye::compute_thermal_energy(
    *params.T_0, debye_temperature, *params.napfu);
  return eta_s * (E_th - E_th_ref) / volume;
}

// Second order (SLB2)
double SLB2::compute_shear_modulus(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return BM2::compute_second_order_shear_modulus(volume, params)
    - compute_shear_modulus_delta(temperature, volume, params);
}

// Third order (SLB3)
double SLB3::compute_shear_modulus(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return BM3::compute_third_order_shear_modulus(volume, params)
    - compute_shear_modulus_delta(temperature, volume, params);
}

double SLB3::compute_molar_heat_capacity_v(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double debye_temperature = compute_debye_temperature(
    *params.V_0 / volume,
    params
  );
  return debye::compute_molar_heat_capacity_v(
    temperature,
    debye_temperature,
    *params.napfu);
}

double SLB3::compute_molar_heat_capacity_p(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double C_v = compute_molar_heat_capacity_v(
    pressure,
    temperature,
    volume,
    params
  );
  double gamma = compute_grueneisen_parameter(
    pressure,
    temperature,
    volume,
    params
  );
  double alpha = compute_thermal_expansivity(
    pressure,
    temperature,
    volume,
    params
  );
  return C_v * (1.0 + gamma * alpha * temperature);
}

double SLB3::compute_thermal_expansivity(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double x = *params.V_0 / volume;
  double debye_temperature = compute_debye_temperature(x, params);
  double C_v = debye::compute_molar_heat_capacity_v(
    temperature,
    debye_temperature,
    *params.napfu
  );
  double gamma_slb = compute_slb_grueneisen_parameter(x, params);
  double K_T = compute_isothermal_bulk_modulus_reuss(
    pressure,
    temperature,
    volume,
    params
  );
  return gamma_slb * C_v / volume / K_T;
}

double SLB3::compute_gibbs_free_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return compute_helmholtz_free_energy(
    pressure,
    temperature,
    volume,
    params
  )
  + pressure * volume;
}

double SLB3::compute_entropy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double debye_temperature = compute_debye_temperature(
    *params.V_0 / volume,
    params
  );
  return debye::compute_entropy(
    temperature,
    debye_temperature,
    *params.napfu);
}

double SLB3::compute_molar_internal_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return compute_helmholtz_free_energy(
    pressure,
    temperature,
    volume,
    params)
  + temperature
  * compute_entropy(
    pressure,
    temperature,
    volume,
    params);
}

double SLB3::compute_helmholtz_free_energy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  double x = *params.V_0 / volume;
  double x_cbrt = std::cbrt(x);
  double x_23 = x_cbrt * x_cbrt;
  double f = 0.5 * (x_23 - 1.0);
  double f2 = f * f;
  double debye_temperature = compute_debye_temperature(x, params);
  double b_iikk = 9.0 * (*params.K_0); // Eq.28
  double b_iikkmm = 27.0 * (*params.K_0) * (*params.Kprime_0 - 4.0); // Eq.29
  double F_quasiharmonic =
    debye::compute_helmholtz_free_energy(
        temperature, debye_temperature, *params.napfu)
    - debye::compute_helmholtz_free_energy(
        *params.T_0, debye_temperature, *params.napfu);
  constexpr double ONE_SIXTH = 1.0 / 6.0;
  return *params.F_0
    + 0.5 * b_iikk * f2 * (*params.V_0)
    + ONE_SIXTH * (*params.V_0) * b_iikkmm * f2 * f
    + F_quasiharmonic;
}

double SLB3::compute_enthalpy(
  double pressure,
  double temperature,
  double volume,
  const MineralParams& params
) const {
  return compute_helmholtz_free_energy(
    pressure,
    temperature,
    volume,
    params)
  + temperature
  * compute_entropy(
    pressure,
    temperature,
    volume,
    params)
  + pressure * volume;
}

double SLB3::compute_debye_temperature(
  double x,
  const MineralParams& params
) {
  double gamma_0 = *params.grueneisen_0;
  double x_cbrt = std::cbrt(x);
  double x_23 = x_cbrt * x_cbrt;
  double f = 0.5 * (x_23 - 1.0);
  // Eq. 47
  double a1_ii = 6.0 * gamma_0;
  double a2_iikk = -12.0 * gamma_0
    + 36.0 * gamma_0 * gamma_0
    - 18.0 * (*params.q_0) * gamma_0;
  // Eq. 41
  double nu_o_nu0_sq = 1.0 + a1_ii * f + 0.5 * a2_iikk * f * f;
  if (nu_o_nu0_sq < 0.0) {
    throw std::logic_error("Volume outside valid range of SLB EoS!");
  }
  return *params.debye_0 * std::sqrt(nu_o_nu0_sq);
}

double SLB3::compute_volume_dependent_q(
  double x,
  const MineralParams& params
) {
  double gamma_0 = *params.grueneisen_0;
  double x_cbrt = std::cbrt(x);
  double x_23 = x_cbrt * x_cbrt;
  double f = 0.5 * (x_23 - 1.0);
  // Eq. 47
  double a1_ii = 6.0 * gamma_0;
  double a2_iikk = -12.0 * gamma_0
    + 36.0 * gamma_0 * gamma_0
    - 18.0 * (*params.q_0) * gamma_0;
  // Eq. 41
  double nu_o_nu0_sq = 1.0 + a1_ii * f + 0.5 * a2_iikk * f * f;
  double two_f_plus1 = 2.0 * f + 1.0;
  constexpr double ONE_SIXTH = 1.0 / 6.0;
  double constexpr ONE_NINTH = 1.0 / 9.0;
  double gr = ONE_SIXTH / nu_o_nu0_sq * two_f_plus1 * (a1_ii + a2_iikk * f);
  double q;
  if (std::abs(gamma_0) < 1.0e-10) {
    q = ONE_NINTH * (18.0 * gr - 6.0);
  } else {
    q = ONE_NINTH * (
      18.0 * gr
      - 6.0
      - 0.5 / nu_o_nu0_sq
      * two_f_plus1 * two_f_plus1
      * a2_iikk
      / gr
    );
  }
  return q;
}

double SLB3::compute_isotropic_eta_s(
  double x,
  const MineralParams& params
) {
  double gamma_0 = *params.grueneisen_0;
  double x_cbrt = std::cbrt(x);
  double x_23 = x_cbrt * x_cbrt;
  double f = 0.5 * (x_23 - 1.0);
  // Eq. 47
  double a2_s = -2.0 * gamma_0 - 2.0 * (*params.eta_s_0);
  double a1_ii = 6.0 * gamma_0;
  double a2_iikk = -12.0 * gamma_0
    + 36.0 * gamma_0 * gamma_0
    - 18.0 * (*params.q_0) * gamma_0;
  // Eq. 41
  double nu_o_nu0_sq = 1.0 + a1_ii * f + 0.5 * a2_iikk * f * f;
  double two_f_plus1 = 2.0 * f + 1.0;
  constexpr double ONE_SIXTH = 1.0 / 6.0;
  double gr = ONE_SIXTH / nu_o_nu0_sq * two_f_plus1 * (a1_ii + a2_iikk * f);
  // Eq. 46 (type in Stixrude 2005)
  return -gr - (0.5 * (1.0 / nu_o_nu0_sq) * two_f_plus1 * two_f_plus1 * a2_s);
}