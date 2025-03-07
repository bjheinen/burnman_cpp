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
#include "burnman/eos/einstein.hpp"
#include "burnman/utils/constants.hpp"

double einstein::compute_thermal_energy_impl(
  double temperature,
  double einstein_temperature,
  double napfu) {
  if (temperature <= constants::precision::double_eps) {
    return 0.0;
  }
  double x = einstein_temperature / temperature;
  // x > ~10 is unlikely, so no need to worry about overflows in expm1(x)
  double exp_term = 1.0 / std::expm1(x);
  // E_th
  return 3.0 * napfu
    * constants::physics::gas_constant
    * einstein_temperature * exp_term;
}

double einstein::compute_molar_heat_capacity_v_impl(
  double temperature,
  double einstein_temperature,
  double napfu) {
  if (temperature < constants::precision::double_eps) {
    return 0.0;
  }
  double x = einstein_temperature / temperature;
  double ex = std::exp(x);
  double exm1 = std::expm1(x);
  // C_v
  return 3.0 * napfu * constants::physics::gas_constant
    * (x * x * ex / (exm1 * exm1));
}

double einstein::compute_helmholtz_free_energy_impl(
  double temperature,
  double einstein_temperature,
  double napfu) {
  double E = compute_thermal_energy_impl(
    temperature, einstein_temperature, napfu);
  double S = compute_entropy_impl(temperature, einstein_temperature, napfu);
  // F
  return E - temperature * S;
}

double einstein::compute_entropy_impl(
  double temperature,
  double einstein_temperature,
  double napfu) {
  if (temperature <= constants::precision::double_eps) {
    return 0.0;
  }
  double x = einstein_temperature / temperature;
  // -x * exp(-x) / (exp(-x) - 1.0) --> x / expm1(x)
  // log(1.0 - exp(-x)) --> log(expm1(x)) - x
  double exm1 = std::expm1(x);
  // S
  return 3.0 * napfu * constants::physics::gas_constant
    * (x / exm1) - (std::log(exm1) - x);
}

double einstein::compute_dmolar_heat_capacity_v_dT_impl(
  double temperature,
  double einstein_temperature,
  double napfu) {
  if (temperature <= constants::precision::double_eps) {
    return 0.0;
  }
  double x = einstein_temperature / temperature;
  // dCvdT
  double ex = std::exp(x);
  double exm1 = std::expm1(x);
  return 3.0 * napfu * constants::physics::gas_constant
    * x * x * ex
    * ((x - 2.0) * ex + (x + 2.0))
    / (temperature * exm1 * exm1 * exm1);
}

double einstein::compute_thermal_energy(
  double temperature,
  double einstein_temperature,
  int napfu) {
  return compute_thermal_energy_impl(
    temperature,
    einstein_temperature,
    static_cast<double>(napfu)
  );
}

double einstein::compute_thermal_energy(
  double temperature,
  double einstein_temperature,
  einstein::ExplicitDouble f) {
  return compute_thermal_energy_impl(
    temperature,
    einstein_temperature,
    f.value
  );
}

double einstein::compute_molar_heat_capacity_v(
  double temperature,
  double einstein_temperature,
  int napfu) {
  return compute_molar_heat_capacity_v_impl(
    temperature,
    einstein_temperature,
    static_cast<double>(napfu)
  );
}

double einstein::compute_molar_heat_capacity_v(
  double temperature,
  double einstein_temperature,
  einstein::ExplicitDouble f) {
  return compute_molar_heat_capacity_v_impl(
    temperature,
    einstein_temperature,
    f.value
  );
}

double einstein::compute_helmholtz_free_energy(
  double temperature,
  double einstein_temperature,
  int napfu) {
  return compute_helmholtz_free_energy_impl(
    temperature,
    einstein_temperature,
    static_cast<double>(napfu)
  );
}

double einstein::compute_helmholtz_free_energy(
  double temperature,
  double einstein_temperature,
  einstein::ExplicitDouble f) {
  return compute_helmholtz_free_energy_impl(
    temperature,
    einstein_temperature,
    f.value
  );
}

double einstein::compute_entropy(
  double temperature,
  double einstein_temperature,
  int napfu) {
  return compute_entropy_impl(
    temperature,
    einstein_temperature,
    static_cast<double>(napfu)
  );
}

double einstein::compute_entropy(
  double temperature,
  double einstein_temperature,
  einstein::ExplicitDouble f) {
  return compute_entropy_impl(
    temperature,
    einstein_temperature,
    f.value
  );
}

double einstein::compute_dmolar_heat_capacity_v_dT(
  double temperature,
  double einstein_temperature,
  int napfu) {
  return compute_dmolar_heat_capacity_v_dT_impl(
    temperature,
    einstein_temperature,
    static_cast<double>(napfu)
  );
}

double einstein::compute_dmolar_heat_capacity_v_dT(
  double temperature,
  double einstein_temperature,
  einstein::ExplicitDouble f) {
  return compute_dmolar_heat_capacity_v_dT_impl(
    temperature,
    einstein_temperature,
    f.value
  );
}
