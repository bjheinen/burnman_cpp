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
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "tolerances.hpp"
#include "burnman/core/mineral.hpp"
#include "burnman/core/equation_of_state.hpp"
#include "burnman/eos/birch_murnaghan.hpp"
#include "burnman/eos/mie_grueneisen_debye.hpp"
#include "burnman/utils/eos.hpp"
#include "burnman/utils/types.hpp"
#include <memory>
#include <typeinfo>
using namespace Catch::Matchers;

TEST_CASE("Set method", "[core][mineral]") {
  Mineral test_mineral;
  test_mineral.params.equation_of_state = EOSType::BM3;
  // Set method should call validate parameters
  REQUIRE_THROWS(test_mineral.set_method(EOSType::Auto));
  test_mineral.params.V_0 = 11.24e-6;
  test_mineral.params.K_0 = 161.0e9;
  test_mineral.params.Kprime_0 = 3.8;
  test_mineral.params.molar_mass = 0.0403;
  test_mineral.params.napfu = 2;
  test_mineral.params.debye_0 = 773.0;
  test_mineral.params.grueneisen_0 = 1.5;
  test_mineral.params.q_0 = 1.5;
  // Test Auto
  REQUIRE_NOTHROW(test_mineral.set_method(EOSType::Auto));
  // Check parameter has not been changed
  REQUIRE(test_mineral.params.equation_of_state == EOSType::BM3);
  // Check eos_method created
  REQUIRE(typeid(*test_mineral.eos_method) == typeid(BM3));
  // Test update EOS
  REQUIRE_NOTHROW(test_mineral.set_method(EOSType::MGD3));
  REQUIRE(test_mineral.params.equation_of_state == EOSType::MGD3);
  REQUIRE(typeid(*test_mineral.eos_method) == typeid(MGD3));
  // Test setting custom EOS (derived from EquationOfState)
  // Using non-derived class will not compile - so no test needed
  class CustomEOS : public EquationOfState{
   public:
    // Helper functions
    bool validate_parameters(
      MineralParams& params [[maybe_unused]]) override {
      return 1;
    }
  };
  CustomEOS custom_eos;
  REQUIRE_NOTHROW(test_mineral.set_method(std::make_shared<CustomEOS>()));
  REQUIRE(test_mineral.params.equation_of_state == EOSType::Custom);
  REQUIRE(typeid(*test_mineral.eos_method) == typeid(CustomEOS));
}

TEST_CASE("Set state", "[core][mineral]") {
  // Set-up mineral
  Mineral test_mineral;
  test_mineral.params.equation_of_state = EOSType::BM3;
  test_mineral.params.V_0 = 11.24e-6;
  test_mineral.params.K_0 = 161.0e9;
  test_mineral.params.Kprime_0 = 3.8;
  test_mineral.params.molar_mass = 0.0403;
  test_mineral.params.napfu = 2;
  test_mineral.params.debye_0 = 773.0;
  test_mineral.params.grueneisen_0 = 1.5;
  test_mineral.params.q_0 = 1.5;
  test_mineral.set_method(EOSType::Auto);
  double test_P = 24e9;
  double test_T = 2000.0;
  // Ensure P/T set correctly
  test_mineral.set_state(test_P, test_T);
  REQUIRE(test_mineral.get_pressure() == test_P);
  REQUIRE(test_mineral.get_temperature() == test_T);
  // Check modifiers still 0
  excesses::Excesses zero_excess;
  excesses::Excesses test_excess = test_mineral.get_property_modifiers();
  CHECK(test_excess.G == zero_excess.G);
  CHECK(test_excess.dGdT == zero_excess.dGdT);
  CHECK(test_excess.dGdP == zero_excess.dGdP);
  CHECK(test_excess.d2GdT2 == zero_excess.d2GdT2);
  CHECK(test_excess.d2GdP2 == zero_excess.d2GdP2);
  CHECK(test_excess.d2GdPdT == zero_excess.d2GdPdT);
  // Define property modifiers and re-check to see if they are computed
  excesses::ExcessParamVector excess_params = {
    excesses::LandauParams{800.0, 1.0e-7, 5.0},
    excesses::LandauSLB2022Params {800.0, 1.0e-8, 5.0},
    excesses::LandauHPParams{298.15, 1.0e-5, 800.0, 1.0e-7, 5.0},
    excesses::LinearParams{1.0e-7, 5.0, 1200.0},
    excesses::BraggWilliamsParams{1, 0.8, 1000.0, 1.0e-7, 1000.0, 1.0e-7},
    excesses::MagneticChsParams{0.4, 800.0, 1.0e-8, 2.2, 1.0e-10},
    excesses::DebyeParams{1.0, 1200.0},
    excesses::DebyeDeltaParams{1.0, 1200.0},
    excesses::EinsteinParams{1.0, 1200.0},
    excesses::EinsteinDeltaParams{1.0, 1200.0}
  };
  REQUIRE_NOTHROW(test_mineral.set_property_modifier_params(excess_params));
  test_mineral.set_state(test_P, test_T);
  CHECK_THAT(test_excess.G,
    !WithinRel(zero_excess.G, tol_rel) &&
    !WithinAbs(zero_excess.G, tol_abs));
  CHECK_THAT(test_excess.dGdT,
    !WithinRel(zero_excess.dGdT, tol_rel) &&
    !WithinAbs(zero_excess.dGdT, tol_abs));
  CHECK_THAT(test_excess.dGdP,
    !WithinRel(zero_excess.dGdP, tol_rel) &&
    !WithinAbs(zero_excess.dGdP, tol_abs));
  CHECK_THAT(test_excess.d2GdT2,
    !WithinRel(zero_excess.d2GdT2, tol_rel) &&
    !WithinAbs(zero_excess.d2GdT2, tol_abs));
  CHECK_THAT(test_excess.d2GdP2,
    !WithinRel(zero_excess.d2GdP2, tol_rel) &&
    !WithinAbs(zero_excess.d2GdP2, tol_abs));
  CHECK_THAT(test_excess.d2GdPdT,
    !WithinRel(zero_excess.d2GdPdT, tol_rel) &&
    !WithinAbs(zero_excess.d2GdPdT, tol_abs));
}

TEST_CASE("Check exceptions", "[core][mineral]") {
  Mineral test_mineral;
  REQUIRE_THROWS(test_mineral.get_molar_mass());
  REQUIRE_THROWS(test_mineral.set_method(EOSType::Auto));
}

TEST_CASE("Check formula", "[core][mineral]") {
  Mineral test_mineral;
  REQUIRE_THROWS(test_mineral.get_formula());
  FormulaMap fm = {
    {"Al", 2.0},
    {"Si", 1.0},
    {"O", 5.0}
  };
  test_mineral.params.formula = fm;
  REQUIRE_NOTHROW(test_mineral.get_formula());
  // TODO (C++20 has ==, but watch doubles -- impement == for FormulaMap)
}

TEST_CASE("Check get/set name", "[core][mineral]") {
  Mineral test_mineral;
  std::string n = "My Mineral!";
  REQUIRE_FALSE(test_mineral.get_name() == n);
  test_mineral.params.name = n;
  REQUIRE(test_mineral.get_name() == n);
}

TEST_CASE("Check py reference values", "[core][mineral]") {
  // Set-up mineral
  Mineral test_mineral;
  test_mineral.params.equation_of_state = EOSType::MGD3;
  test_mineral.params.V_0 = 11.24e-6;
  test_mineral.params.K_0 = 161.0e9;
  test_mineral.params.Kprime_0 = 3.8;
  test_mineral.params.G_0 = 131.0e9;
  test_mineral.params.Gprime_0 = 2.1;
  test_mineral.params.molar_mass = 0.0403;
  test_mineral.params.napfu = 2;
  test_mineral.params.debye_0 = 773.0;
  test_mineral.params.grueneisen_0 = 1.5;
  test_mineral.params.q_0 = 1.5;
  // Set P & T
  double P = 55.0e9;
  double T = 1000.0;

  SECTION("No excess") {
    test_mineral.set_method(EOSType::Auto);
    test_mineral.set_state(P, T);
    // Define reference values
    double ref_Vo = 9.081760466382554e-06;
    double ref_V = 9.081760466382554e-06;
    double ref_m = 0.0403;
    double ref_rho = 4437.465637766628;
    double ref_E = 81891.74478872975;
    double ref_G = 514418.0184366068;
    double ref_F = 14921.192785566382;
    double ref_S = 66.97055200316336;
    double ref_H = 581388.5704397701;
    double ref_KT = 335432555918.34296;
    double ref_KS = 341627283344.15936;
    double ref_invKT = 2.9812252339735266e-12;
    double ref_invKS = 2.9271666777052704e-12;
    double ref_mu = 212373867679.2852;
    double ref_gamma = 1.0894237950619814;
    double ref_alpha = 1.6951968373686158e-05;
    double ref_Cv = 47.40220358077695;
    double ref_Cp = 48.27762168013835;
    double ref_grad = 3.188925030804647e-09;
    double ref_vp = 11865.891735926718;
    double ref_vphi = 8774.225105744063;
    double ref_vs = 6918.039488312308;
    double test_Vo = test_mineral.get_molar_volume_unmodified();
    double test_V = test_mineral.get_molar_volume();
    double test_m = test_mineral.get_molar_mass();
    double test_rho = test_mineral.get_density();
    double test_E = test_mineral.get_molar_internal_energy();
    double test_G = test_mineral.get_molar_gibbs();
    double test_F = test_mineral.get_molar_helmholtz();
    double test_S = test_mineral.get_molar_entropy();
    double test_H = test_mineral.get_molar_enthalpy();
    double test_KT = test_mineral.get_isothermal_bulk_modulus_reuss();
    double test_KS = test_mineral.get_isentropic_bulk_modulus_reuss();
    double test_invKT = test_mineral.get_isothermal_compressibility_reuss();
    double test_invKS = test_mineral.get_isentropic_compressibility_reuss();
    double test_mu = test_mineral.get_shear_modulus();
    double test_gamma = test_mineral.get_grueneisen_parameter();
    double test_alpha = test_mineral.get_thermal_expansivity();
    double test_Cv = test_mineral.get_molar_heat_capacity_v();
    double test_Cp = test_mineral.get_molar_heat_capacity_p();
    double test_grad = test_mineral.get_isentropic_thermal_gradient();
    double test_vp = test_mineral.get_p_wave_velocity();
    double test_vphi = test_mineral.get_bulk_sound_velocity();
    double test_vs = test_mineral.get_shear_wave_velocity();
    CHECK_THAT(test_Vo,
      WithinRel(ref_Vo, tol_rel) || WithinAbs(ref_Vo, tol_abs));
    CHECK_THAT(test_V,
      WithinRel(ref_V, tol_rel) || WithinAbs(ref_V, tol_abs));
    CHECK_THAT(test_m,
      WithinRel(ref_m, tol_rel) || WithinAbs(ref_m, tol_abs));
    CHECK_THAT(test_rho,
      WithinRel(ref_rho, tol_rel) || WithinAbs(ref_rho, tol_abs));
    CHECK_THAT(test_E,
      WithinRel(ref_E, tol_rel) || WithinAbs(ref_E, tol_abs));
    CHECK_THAT(test_G,
      WithinRel(ref_G, tol_rel) || WithinAbs(ref_G, tol_abs));
    CHECK_THAT(test_F,
      WithinRel(ref_F, tol_rel) || WithinAbs(ref_F, tol_abs));
    CHECK_THAT(test_S,
      WithinRel(ref_S, tol_rel) || WithinAbs(ref_S, tol_abs));
    CHECK_THAT(test_H,
      WithinRel(ref_H, tol_rel) || WithinAbs(ref_H, tol_abs));
    CHECK_THAT(test_KT,
      WithinRel(ref_KT, tol_rel) || WithinAbs(ref_KT, tol_abs));
    CHECK_THAT(test_KS,
      WithinRel(ref_KS, tol_rel) || WithinAbs(ref_KS, tol_abs));
    CHECK_THAT(test_invKT,
      WithinRel(ref_invKT, tol_rel) || WithinAbs(ref_invKT, tol_abs));
    CHECK_THAT(test_invKS,
      WithinRel(ref_invKS, tol_rel) || WithinAbs(ref_invKS, tol_abs));
    CHECK_THAT(test_mu,
      WithinRel(ref_mu, tol_rel) || WithinAbs(ref_mu, tol_abs));
    CHECK_THAT(test_gamma,
      WithinRel(ref_gamma, tol_rel) || WithinAbs(ref_gamma, tol_abs));
    CHECK_THAT(test_alpha,
      WithinRel(ref_alpha, tol_rel) || WithinAbs(ref_alpha, tol_abs));
    CHECK_THAT(test_Cv,
      WithinRel(ref_Cv, tol_rel) || WithinAbs(ref_Cv, tol_abs));
    CHECK_THAT(test_Cp,
      WithinRel(ref_Cp, tol_rel) || WithinAbs(ref_Cp, tol_abs));
    CHECK_THAT(test_grad,
      WithinRel(ref_grad, tol_rel) || WithinAbs(ref_grad, tol_abs));
    CHECK_THAT(test_vp,
      WithinRel(ref_vp, tol_rel) || WithinAbs(ref_vp, tol_abs));
    CHECK_THAT(test_vphi,
      WithinRel(ref_vphi, tol_rel) || WithinAbs(ref_vphi, tol_abs));
    CHECK_THAT(test_vs,
      WithinRel(ref_vs, tol_rel) || WithinAbs(ref_vs, tol_abs));
  }

  SECTION("With excess") {
    excesses::ExcessParamVector excess_params = {
      excesses::MagneticChsParams{0.4, 800.0, 1.0e-8, 2.2, 1.0e-10}
    };
    test_mineral.set_property_modifier_params(excess_params);
    test_mineral.set_method(EOSType::Auto);
    test_mineral.set_state(P, T);
    // Define reference values
    double ref_Vo = 9.081760466382554e-06;
    double ref_V = 8.917144831176037e-06;
    double ref_m = 0.0403;
    double ref_rho = 4519.3838120811415;
    double ref_E = 72397.12933495894;
    double ref_G = 509295.0567600913;
    double ref_F = 18852.091045409266;
    double ref_S = 53.54503828954968;
    double ref_H = 562840.095049641;
    double ref_KT = 307455565306.87885;
    double ref_KS = 327195790142.47736;
    double ref_invKT = 1.0 / ref_KT; //3.25250251691452e-12
    double ref_invKS = 1.0 / ref_KS; //3.0562740418039917e-12
    double ref_mu = 212373867679.2852;
    double ref_gamma = 1.7310765377441932;
    double ref_alpha = 3.70897124201741e-05;
    double ref_Cv = 58.741546346304965;
    double ref_Cp = 62.51305502239413;
    double ref_grad = 5.2906442866835055e-09;
    double ref_vp = 11621.274412079585;
    double ref_vphi = 8508.720164047194;
    double ref_vs = 6855.0547115317995;
    double test_Vo = test_mineral.get_molar_volume_unmodified();
    double test_V = test_mineral.get_molar_volume();
    double test_m = test_mineral.get_molar_mass();
    double test_rho = test_mineral.get_density();
    double test_E = test_mineral.get_molar_internal_energy();
    double test_G = test_mineral.get_molar_gibbs();
    double test_F = test_mineral.get_molar_helmholtz();
    double test_S = test_mineral.get_molar_entropy();
    double test_H = test_mineral.get_molar_enthalpy();
    double test_KT = test_mineral.get_isothermal_bulk_modulus_reuss();
    double test_KS = test_mineral.get_isentropic_bulk_modulus_reuss();
    double test_invKT = test_mineral.get_isothermal_compressibility_reuss();
    double test_invKS = test_mineral.get_isentropic_compressibility_reuss();
    double test_mu = test_mineral.get_shear_modulus();
    double test_gamma = test_mineral.get_grueneisen_parameter();
    double test_alpha = test_mineral.get_thermal_expansivity();
    double test_Cv = test_mineral.get_molar_heat_capacity_v();
    double test_Cp = test_mineral.get_molar_heat_capacity_p();
    double test_grad = test_mineral.get_isentropic_thermal_gradient();
    double test_vp = test_mineral.get_p_wave_velocity();
    double test_vphi = test_mineral.get_bulk_sound_velocity();
    double test_vs = test_mineral.get_shear_wave_velocity();
    CHECK_THAT(test_Vo,
      WithinRel(ref_Vo, tol_rel) || WithinAbs(ref_Vo, tol_abs));
    CHECK_THAT(test_V,
      WithinRel(ref_V, tol_rel) || WithinAbs(ref_V, tol_abs));
    CHECK_THAT(test_m,
      WithinRel(ref_m, tol_rel) || WithinAbs(ref_m, tol_abs));
    CHECK_THAT(test_rho,
      WithinRel(ref_rho, tol_rel) || WithinAbs(ref_rho, tol_abs));
    CHECK_THAT(test_E,
      WithinRel(ref_E, tol_rel) || WithinAbs(ref_E, tol_abs));
    CHECK_THAT(test_G,
      WithinRel(ref_G, tol_rel) || WithinAbs(ref_G, tol_abs));
    CHECK_THAT(test_F,
      WithinRel(ref_F, tol_rel) || WithinAbs(ref_F, tol_abs));
    CHECK_THAT(test_S,
      WithinRel(ref_S, tol_rel) || WithinAbs(ref_S, tol_abs));
    CHECK_THAT(test_H,
      WithinRel(ref_H, tol_rel) || WithinAbs(ref_H, tol_abs));
    CHECK_THAT(test_KT,
      WithinRel(ref_KT, tol_rel) || WithinAbs(ref_KT, tol_abs));
    CHECK_THAT(test_KS,
      WithinRel(ref_KS, tol_rel) || WithinAbs(ref_KS, tol_abs));
    CHECK_THAT(test_invKT,
      WithinRel(ref_invKT, tol_rel) || WithinAbs(ref_invKT, tol_abs));
    CHECK_THAT(test_invKS,
      WithinRel(ref_invKS, tol_rel) || WithinAbs(ref_invKS, tol_abs));
    CHECK_THAT(test_mu,
      WithinRel(ref_mu, tol_rel) || WithinAbs(ref_mu, tol_abs));
    CHECK_THAT(test_gamma,
      WithinRel(ref_gamma, tol_rel) || WithinAbs(ref_gamma, tol_abs));
    CHECK_THAT(test_alpha,
      WithinRel(ref_alpha, tol_rel) || WithinAbs(ref_alpha, tol_abs));
    CHECK_THAT(test_Cv,
      WithinRel(ref_Cv, tol_rel) || WithinAbs(ref_Cv, tol_abs));
    CHECK_THAT(test_Cp,
      WithinRel(ref_Cp, tol_rel) || WithinAbs(ref_Cp, tol_abs));
    CHECK_THAT(test_grad,
      WithinRel(ref_grad, tol_rel) || WithinAbs(ref_grad, tol_abs));
    CHECK_THAT(test_vp,
      WithinRel(ref_vp, tol_rel) || WithinAbs(ref_vp, tol_abs));
    CHECK_THAT(test_vphi,
      WithinRel(ref_vphi, tol_rel) || WithinAbs(ref_vphi, tol_abs));
    CHECK_THAT(test_vs,
      WithinRel(ref_vs, tol_rel) || WithinAbs(ref_vs, tol_abs));
  }
}
