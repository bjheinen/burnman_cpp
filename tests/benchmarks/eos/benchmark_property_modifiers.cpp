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
#include "burnman/eos/components/property_modifiers.hpp"
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/eos/components/excess_params.hpp"
#include "burnman/utils/types/mineral_params.hpp"
#include "burnman/core/mineral.hpp"

using namespace burnman;

TEST_CASE("Property modifier benchmarks", "[eos][prop_mod][!benchmark]") {
  // Set-up params for tests
  double P = 1.0e11;
  double P_low = 5.0e10;
  double T = 1000.0;
  double T_high = 2000.0;
  eos::excesses::LandauParams landau_params = {1200.0, 1.0e-7, 5.0};
  eos::excesses::LandauHPParams landau_hp_params = {298.15, 1.0e5, 1200.0, 1.0e-7, 5.0};
  eos::excesses::LandauSLB2022Params landau_slb_params = {1200.0, 1.0e-7, 5.0};
  eos::excesses::LinearParams linear_params = {1.0e-7, 5.0, 1200.0};
  eos::excesses::BraggWilliamsParams bw_params = {1, 0.8, 1000.0, 1.0e-7, 1000.0, 1.0e-7};
  eos::excesses::MagneticChsParams mag_params = {0.4, 1200.0, 1.0e-8, 2.2, 1.0e-10};
  eos::excesses::DebyeParams deb_params = {1.0, 1200.0};
  eos::excesses::DebyeDeltaParams deb_d_params = {1.0, 1200.0};
  eos::excesses::EinsteinParams ein_params = {1.0, 1200.0};
  eos::excesses::EinsteinDeltaParams ein_d_params = {1.0, 1200.0};
  // Set up combined modifiers for integrated test
  eos::excesses::ExcessParamVector excess_params = {
    landau_params,
    landau_hp_params,
    landau_slb_params,
    linear_params,
    bw_params,
    mag_params,
    deb_params,
    deb_d_params,
    ein_params,
    ein_d_params
  };
  // Set up mineral object for integrated test
  Mineral test_mineral;
  test_mineral.params.equation_of_state = types::EOSType::BM3;
  test_mineral.params.V_0 = 11.24e-6;
  test_mineral.params.K_0 = 161.0e9;
  test_mineral.params.Kprime_0 = 3.8;
  test_mineral.set_method(types::EOSType::Auto);
  test_mineral.set_state(P, T);
  test_mineral.set_property_modifier_params(excess_params);

  // Benchmarks for individual modifiers
  BENCHMARK("Landau") {
    return eos::excesses::compute_excesses(P, T, landau_params);
  };

  BENCHMARK("Landau HP") {
    return eos::excesses::compute_excesses(P, T, landau_hp_params);
  };

  BENCHMARK("Landau SLB") {
    return eos::excesses::compute_excesses(P, T, landau_slb_params);
  };

  BENCHMARK("Linear") {
    return eos::excesses::compute_excesses(P, T, linear_params);
  };

  BENCHMARK("Bragg-Williams") {
    return eos::excesses::compute_excesses(P, T, bw_params);
  };

  BENCHMARK("Bragg-Williams (Low P - High T)") {
    return eos::excesses::compute_excesses(P_low, T_high, bw_params);
  };

  BENCHMARK("Magnetic Chs") {
    return eos::excesses::compute_excesses(P, T, mag_params);
  };

  BENCHMARK("Debye") {
    return eos::excesses::compute_excesses(P, T, deb_params);
  };

  BENCHMARK("Debye delta") {
    return eos::excesses::compute_excesses(P, T, deb_d_params);
  };

  BENCHMARK("Einstein") {
    return eos::excesses::compute_excesses(P, T, ein_params);
  };

  BENCHMARK("Einstein delta") {
    return eos::excesses::compute_excesses(P, T, ein_d_params);
  };

  // Integrated benchmark calling all modifiers from Mineral
  BENCHMARK("Combined modifiers") {
    test_mineral.set_state(P, T);
    return test_mineral.get_property_modifiers();
  };

}
