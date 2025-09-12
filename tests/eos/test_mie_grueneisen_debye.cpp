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
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "burnman/eos/mie_grueneisen_debye.hpp"
#include <cmath>
#include <map>
#include <tuple>
#include "burnman/utils/types/mineral_params.hpp"
#include "tolerances.hpp"

using namespace Catch::Matchers;
using namespace burnman;

TEST_CASE("Test validate parameters", "[mgd][eos]") {
  eos::MGD3 mgd;
  SECTION("Missing V_0") {
    types::MineralParams params;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing K_0") {
    types::MineralParams params;
    params.V_0 = 11.24e-6;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing Kprime_0") {
    types::MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing molar mass") {
    types::MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing napfu") {
    types::MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing debye_0") {
    types::MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing grueneisen_0") {
    types::MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.q_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }
  SECTION("Missing q_0") {
    types::MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    REQUIRE_THROWS(mgd.validate_parameters(params));
  }

  SECTION("Optional Parameters") {
    types::MineralParams params;
    params.V_0 = 11.24e-6;
    params.K_0 = 161.0e9;
    params.Kprime_0 = 3.8;
    params.molar_mass = 0.0403;
    params.napfu = 2;
    params.debye_0 = 773.0;
    params.grueneisen_0 = 1.5;
    params.q_0 = 1.5;
    REQUIRE_FALSE(params.T_0.has_value());
    REQUIRE_FALSE(params.P_0.has_value());
    REQUIRE_FALSE(params.E_0.has_value());
    REQUIRE_FALSE(params.F_0.has_value());
    REQUIRE_FALSE(params.G_0.has_value());
    REQUIRE_FALSE(params.Gprime_0.has_value());
    mgd.validate_parameters(params);
    REQUIRE(*params.T_0 == 300.0);
    REQUIRE(*params.P_0 == 0);
    REQUIRE(*params.E_0 == 0);
    REQUIRE(*params.F_0 == 0);
    REQUIRE(std::isnan(*params.G_0));
    REQUIRE(std::isnan(*params.Gprime_0));
  }
  // TODO: Check warnings
}

TEST_CASE("Check reference conditions", "[mgd][eos]") {
  // Set up test params
  types::MineralParams params;
  params.E_0 = 0.0;
  params.T_0 = 300.0;
  params.P_0 = 0.0;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.G_0 = 131.0e9;
  params.Gprime_0 = 2.1;
  params.debye_0 = 773.0;
  params.grueneisen_0 = 1.5;
  params.q_0 = 1.5;
  params.molar_mass = 0.0403;
  params.napfu = 2;

  eos::MGD2 mgd2;
  eos::MGD3 mgd3;
  REQUIRE_NOTHROW(mgd3.validate_parameters(params));

  SECTION("Independent functions") {
    double T = *params.T_0;
    double P = *params.P_0;
    double V = *params.V_0;
    CHECK_THAT(mgd3.compute_pressure(T, V, params),
      WithinRel(*params.P_0, tol_rel) || WithinAbs(*params.P_0, tol_abs));
    CHECK_THAT(mgd3.compute_volume(P, T, params),
      WithinRel(*params.V_0, tol_rel) || WithinAbs(*params.V_0, tol_abs));
  }
  SECTION("V dependent functions") {
    auto P = GENERATE(0.0, 10.0, 25.e9);
    auto T = GENERATE(300.0, 2000.0);
    double V = *params.V_0;
    CHECK_THAT(mgd3.compute_grueneisen_parameter(P, T, V, params),
      WithinRel(*params.grueneisen_0, tol_rel) ||
      WithinAbs(*params.grueneisen_0, tol_abs));
  }
  SECTION("V-T dependent functions") {
    auto P = GENERATE(0.0, 10.0, 25.e9);
    double T = *params.T_0;
    double V = *params.V_0;
    CHECK_THAT(mgd3.compute_isothermal_bulk_modulus_reuss(P, T, V, params),
      WithinRel(*params.K_0, tol_rel) || WithinAbs(*params.K_0, tol_abs));
    CHECK_THAT(mgd3.compute_shear_modulus(P, T, V, params),
      WithinRel(*params.G_0, tol_rel) || WithinAbs(*params.G_0, tol_abs));
    CHECK_THAT(mgd2.compute_shear_modulus(P, T, V, params),
      WithinRel(*params.G_0, tol_rel) || WithinAbs(*params.G_0, tol_abs));
    CHECK_THAT(mgd3.compute_helmholtz_free_energy(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
  }
  SECTION("P-T-V dependent functions") {
    types::MineralParams params_b = params;
    params_b.T_0 = 500;
    params_b.P_0 = 1.e6;
    double P = *params.P_0;
    double T = *params.T_0;
    double V = *params.V_0;
    double P2 = *params_b.P_0;
    double T2 = *params_b.T_0;
    CHECK_THAT(mgd3.compute_gibbs_free_energy(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
    CHECK_THAT(mgd3.compute_gibbs_free_energy(P2, T2, V, params_b),
      WithinRel(P2*V, tol_rel) || WithinAbs(P2*V, tol_abs));
  }
  SECTION("Zero T known functions") {
    auto P = GENERATE(0.0, 10.0, 25.e9);
    double T = 0.0;
    double V = *params.V_0;
    CHECK_THAT(mgd3.compute_molar_heat_capacity_v(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
    CHECK_THAT(mgd3.compute_molar_heat_capacity_p(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
    CHECK_THAT(mgd3.compute_thermal_expansivity(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
    // At T_ref: K_T = K_0; At T=0: K_S = K_T
    CHECK_THAT(mgd3.compute_isentropic_bulk_modulus_reuss(P, T, V, params),
      !WithinRel(*params.K_0, tol_rel) && !WithinAbs(*params.K_0, tol_abs));
    // Copy params to set T_0 = 0
    types::MineralParams params_T0 = params;
    params_T0.T_0 = 0;
    CHECK_THAT(mgd3.compute_isentropic_bulk_modulus_reuss(P, T, V, params_T0),
      WithinRel(*params.K_0, tol_rel) || WithinAbs(*params.K_0, tol_abs));
  }
}

TEST_CASE("Shear modulus expansion", "[mgd][eos]") {
  // Set up test P, T, V
  double P = 25.e9;
  double T = 2000.0;
  double V = 10.0e-6;
  // Set up test params
  types::MineralParams params;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.debye_0 = 773.0;
  params.grueneisen_0 = 1.5;
  params.q_0 = 1.5;
  params.molar_mass = 0.0403;
  params.napfu = 2;
  eos::MGD2 mgd2;
  eos::MGD3 mgd3;
  // Validating should set G0/G' to nan
  mgd3.validate_parameters(params);
  REQUIRE_NOTHROW(mgd2.compute_shear_modulus(P, T, V, params));
  REQUIRE_NOTHROW(mgd3.compute_shear_modulus(P, T, V, params));
  CHECK(std::isnan(mgd2.compute_shear_modulus(P, T, V, params)));
  CHECK(std::isnan(mgd3.compute_shear_modulus(P, T, V, params)));
  // Set G0/G'0
  params.G_0 = 131.0e9;
  params.Gprime_0 = 2.1;
  CHECK_FALSE(std::isnan(mgd3.compute_shear_modulus(P, T, V, params)));
  // Ensure expansion order different
  CHECK_FALSE(
    mgd2.compute_shear_modulus(P, T, V, params) ==
    mgd3.compute_shear_modulus(P, T, V, params)
  );
}

TEST_CASE("Test MGD volume", "[mgd][eos]") {
  // Set up test params
  types::MineralParams params;
  params.T_0 = 300.0;
  params.P_0 = 0.0;
  params.E_0 = 0.0;
  params.F_0 = 0.0;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.G_0 = 131.0e9;
  params.Gprime_0 = 2.1;
  params.debye_0 = 773.0;
  params.grueneisen_0 = 1.5;
  params.q_0 = 1.4;
  params.molar_mass = 0.0403;
  params.napfu = 2;
  eos::MGD3 mgd;
  double T_a = 800.0;
  double T_b = 2000.0;
  double V_a = 0.9 * (*params.V_0);
  double V_b = 0.5 * (*params.V_0);
  double P_aa = mgd.compute_pressure(T_a, V_a, params);
  double P_ab = mgd.compute_pressure(T_a, V_b, params);
  double P_ba = mgd.compute_pressure(T_b, V_a, params);
  double P_bb = mgd.compute_pressure(T_b, V_b, params);
  CHECK_THAT(mgd.compute_volume(P_aa, T_a, params),
    WithinRel(V_a, tol_rel) || WithinAbs(V_a, tol_abs));
  CHECK_THAT(mgd.compute_volume(P_ab, T_a, params),
    WithinRel(V_b, tol_rel) || WithinAbs(V_b, tol_abs));
  CHECK_THAT(mgd.compute_volume(P_ba, T_b, params),
    WithinRel(V_a, tol_rel) || WithinAbs(V_a, tol_abs));
  CHECK_THAT(mgd.compute_volume(P_bb, T_b, params),
    WithinRel(V_b, tol_rel) || WithinAbs(V_b, tol_abs));
}

TEST_CASE("MGD python reference values", "[mgd][eos]") {
  // Set up test params
  types::MineralParams params;
  params.T_0 = 300.0;
  params.P_0 = 0.0;
  params.E_0 = 0.0;
  params.F_0 = 0.0;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.G_0 = 131.0e9;
  params.Gprime_0 = 2.1;
  params.debye_0 = 773.0;
  params.grueneisen_0 = 1.5;
  params.q_0 = 1.4;
  params.molar_mass = 0.0403;
  params.napfu = 2;

  SECTION("V dependent functions") {
    struct TestData {
      double input;
      double expected_gamma;
    };
    eos::MGD3 mgd3;
    // P & T unused
    double P = 2.0e9;
    double T = 500.0;
    // TODO: Data
    auto test_data = GENERATE(
      TestData{0.99, 1.4790420843377579},
      TestData{0.98, 1.4581686774325902},
      TestData{0.95, 1.3960607156293774},
      TestData{0.80, 1.0975321246255834},
      TestData{0.40, 0.4158869058930879}
    );
    CAPTURE(test_data.input);
    CHECK_THAT(mgd3.compute_grueneisen_parameter(
        P, T, test_data.input*(*params.V_0), params
      ),
      WithinRel(test_data.expected_gamma, tol_rel) ||
      WithinAbs(test_data.expected_gamma, tol_abs));
  }
  SECTION("T-V dependent functions") {
    eos::MGD3 mgd3;
    eos::MGD2 mgd2;
    double P = 2.e9;
    double T1 = 300, T2 = 800, T3 = 2500;
    double x1 = 0.99, x2 = 0.80, x3 = 0.40;
    // Reference data
    std::map<std::tuple<double, double>, double> ref_P = {
      {{T1, x1}, 1649296384.9124484},
      {{T1, x2}, 54834505788.6455},
      {{T1, x3}, 818147973581.2825},
      {{T2, x1}, 4585804419.443684},
      {{T2, x2}, 57322435755.93491},
      {{T2, x3}, 819582614695.9464},
      {{T3, x1}, 15687395638.655682},
      {{T3, x2}, 67406443542.81873},
      {{T3, x3}, 826913879314.8103}
    };
    std::map<std::tuple<double, double>, double> ref_KT = {
      {{T1, x1}, 167236609646.01492},
      {{T1, x2}, 349338640343.3208},
      {{T1, x3}, 2384291228322.9097},
      {{T2, x1}, 165062573508.1445},
      {{T2, x2}, 347331538937.907},
      {{T2, x3}, 2383244538129.3477},
      {{T3, x1}, 160124964007.1976},
      {{T3, x2}, 342727116234.12994},
      {{T3, x3}, 2379913952857.0317}
    };
    std::map<std::tuple<double, double>, double> ref_KS = {
      {{T1, x1}, 169374921707.40085},
      {{T1, x2}, 350523917990.9041},
      {{T1, x3}, 2384466361991.622},
      {{T2, x1}, 172543537457.656},
      {{T2, x2}, 352259329066.90717},
      {{T2, x3}, 2384489154000.001},
      {{T3, x1}, 184522621586.30026},
      {{T3, x2}, 359293248443.23315},
      {{T3, x3}, 2384605625111.078}
    };
    std::map<std::tuple<double, double>, double> ref_G2 = {
      {{T1, x1}, 134440941855.59421},
      {{T1, x2}, 231811619814.15884},
      {{T1, x3}, 1299855849555.416},
      {{T2, x1}, 129612710531.43448},
      {{T2, x2}, 227621843010.16327},
      {{T2, x3}, 1297506266101.6821},
      {{T3, x1}, 113328235367.81192},
      {{T3, x2}, 212758380043.6364},
      {{T3, x3}, 1286710397395.656}
    };
    std::map<std::tuple<double, double>, double> ref_G3 = {
      {{T1, x1}, 134430413937.49939},
      {{T1, x2}, 223263157637.72803},
      {{T1, x3}, 551942177125.8463},
      {{T2, x1}, 129602182613.33966},
      {{T2, x2}, 219073380833.73245},
      {{T2, x3}, 549592593672.1124},
      {{T3, x1}, 113317707449.7171},
      {{T3, x2}, 204209917867.20563},
      {{T3, x3}, 538796724966.08624}
    };
    std::map<std::tuple<double, double>, double> ref_Cv = {
      {{T1, x1}, 36.256867025689054},
      {{T1, x2}, 29.493115813004234},
      {{T1, x3}, 15.17484275132262},
      {{T2, x1}, 47.56724708932584},
      {{T2, x2}, 45.98159509496895},
      {{T2, x3}, 40.44093205026483},
      {{T3, x1}, 49.64191948730763},
      {{T3, x2}, 49.465560088421185},
      {{T3, x3}, 48.782434047774686}
    };
    std::map<std::tuple<double, double>, double> ref_Cp = {
      {{T1, x1}, 36.720452697709064},
      {{T1, x2}, 29.59318356072428},
      {{T1, x3}, 15.17595739111645},
      {{T2, x1}, 49.72308928352981},
      {{T2, x2}, 46.633962142077735},
      {{T2, x3}, 40.46205175705467},
      {{T3, x1}, 57.20567796013561},
      {{T3, x2}, 51.856538127235964},
      {{T3, x3}, 48.87860189116753}
    };
    std::map<std::tuple<double, double>, double> ref_S = {
      {{T1, x1}, 26.1688215204046},
      {{T1, x2}, 17.17207187986112},
      {{T1, x3}, 6.291121348908287},
      {{T2, x1}, 68.66122959604459},
      {{T2, x2}, 55.8970670193532},
      {{T2, x3}, 34.68258890478288},
      {{T3, x1}, 124.44699145969541},
      {{T3, x2}, 110.94109298759564},
      {{T3, x3}, 86.99709822571647}
    };
    std::map<std::tuple<double, double>, double> ref_F = {
      {{T1, x1}, 91.94810310163791},
      {{T1, x2}, 51536.61471078331},
      {{T1, x3}, 1321858.6932650995},
      {{T2, x1}, -24893.517554948095},
      {{T2, x2}, 32354.00995392558},
      {{T2, x3}, 1311509.3346771272},
      {{T3, x1}, -197558.98924149506},
      {{T3, x2}, -117663.53453757522},
      {{T3, x3}, 1201018.258956775}
    };
    std::map<std::tuple<double, double>, double> ref_E = {
      {{T1, x1}, 7942.594559223017},
      {{T1, x2}, 56688.23627474165},
      {{T1, x3}, 1323746.0296697721},
      {{T2, x1}, 30035.46612188758},
      {{T2, x2}, 77071.66356940813},
      {{T2, x3}, 1339255.4058009535},
      {{T3, x1}, 113558.48940774347},
      {{T3, x2}, 159689.1979314139},
      {{T3, x3}, 1418511.0045210663}
    };
    std::map<std::tuple<double, double>, double> ref_alpha = {
      {{T1, x1}, 2.8816281951684226e-05},
      {{T1, x2}, 1.030469099558686e-05},
      {{T1, x3}, 5.887268495167576e-07},
      {{T2, x1}, 3.830349751433756e-05},
      {{T2, x2}, 1.615848929441986e-05},
      {{T2, x3}, 1.5696451926796345e-06},
      {{T3, x1}, 4.120676741826695e-05},
      {{T3, x2}, 1.7616328765565324e-05},
      {{T3, x3}, 1.8960559789244369e-06}
    };
    // Generate combos
    // double T1 = 300, T2 = 800, T3 = 2500;
    // double x1 = 0.99, x2 = 0.80, x3 = 0.40;
    auto T = GENERATE(300.0, 800.0, 2500.0);
    auto x = GENERATE(0.99, 0.8, 0.4);
    CAPTURE(T);
    CAPTURE(x);
    double V = *params.V_0 * x;
    auto key = std::make_tuple(T, x);
    CHECK_THAT(mgd3.compute_pressure(
        T, V, params
      ),
      WithinRel(ref_P[key], tol_rel) ||
      WithinAbs(ref_P[key], tol_abs));
    CHECK_THAT(mgd3.compute_isothermal_bulk_modulus_reuss(
        P, T, V, params
      ),
      WithinRel(ref_KT[key], tol_rel) ||
      WithinAbs(ref_KT[key], tol_abs));
    CHECK_THAT(mgd3.compute_isentropic_bulk_modulus_reuss(
        P, T, V, params
      ),
      WithinRel(ref_KS[key], tol_rel) ||
      WithinAbs(ref_KS[key], tol_abs));
    CHECK_THAT(mgd2.compute_shear_modulus(
        P, T, V, params
      ),
      WithinRel(ref_G2[key], tol_rel) ||
      WithinAbs(ref_G2[key], tol_abs));
    CHECK_THAT(mgd3.compute_shear_modulus(
        P, T, V, params
      ),
      WithinRel(ref_G3[key], tol_rel) ||
      WithinAbs(ref_G3[key], tol_abs));
    CHECK_THAT(mgd3.compute_molar_heat_capacity_v(
        P, T, V, params
      ),
      WithinRel(ref_Cv[key], tol_rel) ||
      WithinAbs(ref_Cv[key], tol_abs));
    CHECK_THAT(mgd3.compute_molar_heat_capacity_p(
        P, T, V, params
      ),
      WithinRel(ref_Cp[key], tol_rel) ||
      WithinAbs(ref_Cp[key], tol_abs));
    CHECK_THAT(mgd3.compute_entropy(
        P, T, V, params
      ),
      WithinRel(ref_S[key], tol_rel) ||
      WithinAbs(ref_S[key], tol_abs));
    CHECK_THAT(mgd3.compute_helmholtz_free_energy(
        P, T, V, params
      ),
      WithinRel(ref_F[key], tol_rel) ||
      WithinAbs(ref_F[key], tol_abs));
    CHECK_THAT(mgd3.compute_molar_internal_energy(
        P, T, V, params
      ),
      WithinRel(ref_E[key], tol_rel) ||
      WithinAbs(ref_E[key], tol_abs));
    CHECK_THAT(mgd3.compute_thermal_expansivity(
        P, T, V, params
      ),
      WithinRel(ref_alpha[key], tol_rel) ||
      WithinAbs(ref_alpha[key], tol_abs));
  }
  SECTION("P-T-V dependent functions") {
    eos::MGD3 mgd3;
    double P1 = 1.e9, P2 = 54.e9;
    double T1 = 800, T2 = 2500;
    double x1 = 0.99, x2 = 0.4;
    std::map<std::tuple<double, double, double>, double> ref_G = {
        {{P1, T1, x1}, -13765.917554948095},
        {{P1, T1, x2}, 1316005.3346771272},
        {{P1, T2, x1}, -186431.38924149505},
        {{P1, T2, x2}, 1205514.258956775},
        {{P2, T1, x1}, 575996.8824450519},
        {{P2, T1, x2}, 1554293.3346771272},
        {{P2, T2, x1}, 403331.410758505},
        {{P2, T2, x2}, 1443802.258956775}
    };
    std::map<std::tuple<double, double, double>, double> ref_H = {
        {{P1, T1, x1}, 41163.06612188758},
        {{P1, T1, x2}, 1343751.4058009535},
        {{P1, T2, x1}, 124686.08940774348},
        {{P1, T2, x2}, 1423007.0045210663},
        {{P2, T1, x1}, 630925.8661218876},
        {{P2, T1, x2}, 1582039.4058009535},
        {{P2, T2, x1}, 714448.8894077435},
        {{P2, T2, x2}, 1661295.0045210663}
    };
    // Generate combos
    auto P = GENERATE(1.e9, 54.e9);
    auto T = GENERATE(800.0, 2500.0);
    auto x = GENERATE(0.99, 0.4);
    CAPTURE(P);
    CAPTURE(T);
    CAPTURE(x);
    double V = *params.V_0 * x;
    auto key = std::make_tuple(P, T, x);
    CHECK_THAT(mgd3.compute_gibbs_free_energy(
        P, T, V, params
      ),
      WithinRel(ref_G[key], tol_rel) ||
      WithinAbs(ref_G[key], tol_abs));
    CHECK_THAT(mgd3.compute_enthalpy(
        P, T, V, params
      ),
      WithinRel(ref_H[key], tol_rel) ||
      WithinAbs(ref_H[key], tol_abs));
  }
}
