/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "burnman/core/mineral.hpp"
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/utils/types/mineral_params.hpp"
#include "eos_names.hpp"

using namespace burnman;

TEST_CASE("Mineral per-property benchmarks", "[core][mineral][!benchmark]") {
  // Generate eos type from list of options
  auto eos_type = GENERATE(
    types::EOSType::MT,
    types::EOSType::Vinet,
    types::EOSType::BM2,
    types::EOSType::BM3,
    types::EOSType::MGD2,
    types::EOSType::MGD3,
    types::EOSType::SLB2,
    types::EOSType::SLB3,
    types::EOSType::SLB3Conductive,
    types::EOSType::HPTMT
  );
  // Capture var and convert to string to add to benchmark names
  CAPTURE(eos_type);
  std::string eos_name = std::string(" [") + eos_string(eos_type) + "]";
  // Set-up mineral object
  Mineral test_mineral;
  // Set mineral EOS from generated type
  test_mineral.params.equation_of_state = eos_type;
  // Set-up common mineral params
  test_mineral.params.T_0 = 300.0;
  test_mineral.params.P_0 = 0.0;
  test_mineral.params.E_0 = 0.0;
  test_mineral.params.F_0 = 0.0;
  test_mineral.params.H_0 = -1443030.0;
  test_mineral.params.V_0 = 11.24e-6;
  test_mineral.params.K_0 = 161.0e9;
  test_mineral.params.Kprime_0 = 3.8;
  test_mineral.params.Kdprime_0 = -1.6e-11;
  test_mineral.params.G_0 = 131.0e9;
  test_mineral.params.Gprime_0 = 2.1;
  test_mineral.params.molar_mass = 0.0403;
  test_mineral.params.napfu = 2;
  test_mineral.params.debye_0 = 773.0;
  test_mineral.params.grueneisen_0 = 1.5;
  test_mineral.params.q_0 = 1.5;
  test_mineral.params.Cp = types::CpParams{149.3, 0.002918, -2983000.0, -799.1};
  test_mineral.params.S_0 = 62.6;
  test_mineral.params.a_0 = 1.87e-05;
  test_mineral.params.eta_s_0 = 2.565;
  test_mineral.params.bel_0 = 0.00411;
  test_mineral.params.gel = 1.476;
  // Set eos from parameter
  test_mineral.set_method(types::EOSType::Auto);
  // Set mineral state
  double P = 55.0e9;
  double T = 1000.0;
  test_mineral.set_state(P, T);
  // Run individual benchmarks

  // Benchmark cache clear as used in every benchmark
  BENCHMARK("reset_cache" + eos_name) {
    test_mineral.reset_cache();
    return 0;
  };

  BENCHMARK("get_molar_volume" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_molar_volume();
  };

  BENCHMARK("get_density" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_density();
  };

  BENCHMARK("get_molar_internal_energy" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_molar_internal_energy();
  };

  BENCHMARK("get_molar_gibbs" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_molar_gibbs();
  };

  BENCHMARK("get_molar_helmholtz" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_molar_helmholtz();
  };

  BENCHMARK("get_molar_entropy" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_molar_entropy();
  };

  BENCHMARK("get_molar_enthalpy" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_molar_enthalpy();
  };

  BENCHMARK("get_isothermal_bulk_modulus_reuss" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_isothermal_bulk_modulus_reuss();
  };

  BENCHMARK("get_isentropic_bulk_modulus_reuss" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_isentropic_bulk_modulus_reuss();
  };

  BENCHMARK("get_isothermal_compressibility_reuss" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_isothermal_compressibility_reuss();
  };

  BENCHMARK("get_isentropic_compressibility_reuss" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_isentropic_compressibility_reuss();
  };

  BENCHMARK("get_shear_modulus" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_shear_modulus();
  };

  BENCHMARK("get_grueneisen_parameter" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_grueneisen_parameter();
  };

  BENCHMARK("get_thermal_expansivity" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_thermal_expansivity();
  };

  BENCHMARK("get_molar_heat_capacity_v" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_molar_heat_capacity_v();
  };

  BENCHMARK("get_molar_heat_capacity_p" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_molar_heat_capacity_p();
  };

  BENCHMARK("get_isentropic_thermal_gradient" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_isentropic_thermal_gradient();
  };

  BENCHMARK("get_p_wave_velocity" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_p_wave_velocity();
  };

  BENCHMARK("get_bulk_sound_velocity" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_bulk_sound_velocity();
  };

  BENCHMARK("get_shear_wave_velocity" + eos_name) {
    test_mineral.reset_cache();
    return test_mineral.get_shear_wave_velocity();
  };
}
