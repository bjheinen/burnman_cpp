/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */

#ifndef TESTS_SOLUTION_FIXTURES_HPP_INCLUDED
#define TESTS_SOLUTION_FIXTURES_HPP_INCLUDED

#include <memory>
#include <Eigen/Dense>
#include "burnman/core/mineral.hpp"
#include "burnman/core/solution_model.hpp"
#include "burnman/core/solution.hpp"
#include "burnman/core/assemblage.hpp"
#include "burnman/utils/types.hpp"

struct BridgmaniteFixture {
  // Declare variables to use
  Mineral mg_si_perovskite;
  Mineral fe_si_perovskite;
  Mineral al_al_perovskite;
  PairedEndmemberList bdg_endmembers;
  std::shared_ptr<IdealSolution> bdg_solution_model;
  Eigen::ArrayXd molar_fractions;
  double P;
  double T;

  BridgmaniteFixture() {
    // Set MgPv params
    mg_si_perovskite.params.name = "MgSiO3 perovskite";
    mg_si_perovskite.params.formula = FormulaMap{
      {"Mg", 1.0},
      {"Si", 1.0},
      {"O", 3.0}
    };
    mg_si_perovskite.params.napfu = 5;
    mg_si_perovskite.params.molar_mass = 0.1003887;
    mg_si_perovskite.params.F_0 = -1368000.0;
    mg_si_perovskite.params.V_0 = 2.4445e-05;
    mg_si_perovskite.params.K_0 = 2.51e11;
    mg_si_perovskite.params.Kprime_0 = 4.1;
    mg_si_perovskite.params.G_0 = 1.73e11;
    mg_si_perovskite.params.Gprime_0 = 1.7;
    mg_si_perovskite.params.debye_0 = 905.0;
    mg_si_perovskite.params.grueneisen_0 = 1.57;
    mg_si_perovskite.params.q_0 = 1.1;
    mg_si_perovskite.params.eta_s_0 = 2.3;
    mg_si_perovskite.params.equation_of_state = EOSType::SLB3;
    // Need to set method/name etc.

    // Set FePv params
    fe_si_perovskite.params.name = "FeSiO3 perovskite";
    fe_si_perovskite.params.formula = FormulaMap{
      {"Fe", 1.0},
      {"Si", 1.0},
      {"O", 3.0}
    };
    fe_si_perovskite.params.napfu = 5;
    fe_si_perovskite.params.molar_mass = 0.1319287;
    fe_si_perovskite.params.F_0 = -1043000.0;
    fe_si_perovskite.params.V_0 = 2.534e-05;
    fe_si_perovskite.params.K_0 = 2.72e11;
    fe_si_perovskite.params.Kprime_0 = 4.1;
    fe_si_perovskite.params.G_0 = 1.33e11;
    fe_si_perovskite.params.Gprime_0 = 1.4;
    fe_si_perovskite.params.debye_0 = 871.0;
    fe_si_perovskite.params.grueneisen_0 = 1.57;
    fe_si_perovskite.params.q_0 = 1.1;
    fe_si_perovskite.params.eta_s_0 = 2.3;
    fe_si_perovskite.params.equation_of_state = EOSType::SLB3;

    // Set AlPv params
    al_al_perovskite.params.name = "AlAlO3 perovskite";
    al_al_perovskite.params.formula = FormulaMap{
      {"Al", 2.0},
      {"O", 3.0}
    };
    al_al_perovskite.params.napfu = 5;
    al_al_perovskite.params.molar_mass = 0.1019612;
    al_al_perovskite.params.F_0 = -1533878.0;
    al_al_perovskite.params.V_0 = 2.494e-05;
    al_al_perovskite.params.K_0 = 2.58e11;
    al_al_perovskite.params.Kprime_0 = 4.1;
    al_al_perovskite.params.G_0 = 1.71e11;
    al_al_perovskite.params.Gprime_0 = 1.5;
    al_al_perovskite.params.debye_0 = 886.0;
    al_al_perovskite.params.grueneisen_0 = 1.57;
    al_al_perovskite.params.q_0 = 1.1;
    al_al_perovskite.params.eta_s_0 = 2.5;
    al_al_perovskite.params.equation_of_state = EOSType::SLB3;

    // Paired list for solution model setup
    bdg_endmembers = {
      {mg_si_perovskite, "[Mg][Si]O3"},
      {fe_si_perovskite, "[Fe][Si]O3"},
      {al_al_perovskite, "[Al][Al]O3"},
    };
    // Make solution model
    bdg_solution_model = std::make_shared<IdealSolution>(bdg_endmembers);

    // Store molar fractions here
    molar_fractions.resize(3);
    molar_fractions << 0.88, 0.07, 0.05;

    // Store P, T
    P = 40.e9;
    T= 2000.0;
  }
};

struct FerropericlaseFixture {
  // Declare variables to use
  Mineral periclase;
  Mineral wuestite;
  excesses::ExcessParamVector excess_params_wuestite;
  PairedEndmemberList fp_endmembers;
  std::shared_ptr<SymmetricRegularSolution> fp_solution_model;
  Eigen::ArrayXd molar_fractions;

  FerropericlaseFixture() {
    // Periclase
    periclase.params.name = "Periclase";
    periclase.params.formula = FormulaMap{
          {"Mg", 1.0},
          {"O", 1.0}
        };
    periclase.params.napfu = 2;
    periclase.params.molar_mass = 0.0403044;
    periclase.params.F_0 = -569444.6;
    periclase.params.V_0 = 1.1244e-05;
    periclase.params.K_0 = 1.613836e11;
    periclase.params.Kprime_0 = 3.84045;
    periclase.params.G_0 = 1.309e11;
    periclase.params.Gprime_0 = 2.1438;
    periclase.params.debye_0 = 767.0977;
    periclase.params.grueneisen_0 = 1.36127;
    periclase.params.q_0 = 1.7217;
    periclase.params.eta_s_0 = 2.81765;
    periclase.params.equation_of_state = EOSType::SLB3;
    // Wuestite
    wuestite.params.name = "Wuestite";
    wuestite.params.formula = FormulaMap{
          {"Fe", 1.0},
          {"O", 1.0}
        };
    wuestite.params.napfu = 2;
    wuestite.params.molar_mass = 0.0718444;
    wuestite.params.F_0 = -242146.0;
    wuestite.params.V_0 = 1.2264e-05;
    wuestite.params.K_0 = 1.794442e11;
    wuestite.params.Kprime_0 = 4.9376;
    wuestite.params.G_0 = 59000000000.0;
    wuestite.params.Gprime_0 = 1.44673;
    wuestite.params.debye_0 = 454.1592;
    wuestite.params.grueneisen_0 = 1.53047;
    wuestite.params.q_0 = 1.7217;
    wuestite.params.eta_s_0 = -0.05731;
    wuestite.params.equation_of_state = EOSType::SLB3;
    excess_params_wuestite = { excesses::LinearParams{1.0, 2.0, 3.0} };
    wuestite.set_property_modifier_params(excess_params_wuestite);

    // Paired list for solution model setup
    fp_endmembers = {
      {periclase, "[Mg]O"},
      {wuestite, "[Fe]O"}
    };
    // Make solution model
    fp_solution_model = std::make_shared<SymmetricRegularSolution>(
      fp_endmembers,
      std::vector<std::vector<double>>{{13.0e3}}
    );
    // Store molar fractions here
    molar_fractions.resize(2);
    molar_fractions << 0.9, 0.1;
  }
};

struct CaPerovskiteFixture {
  // Declare variables to use
  Mineral ca_perovskite;
  CaPerovskiteFixture() {
    ca_perovskite.params.name = "Ca-perovskite";
    ca_perovskite.params.formula = FormulaMap{
          {"Ca", 1.0},
          {"Si", 1.0},
          {"O", 3.0}
        };
    ca_perovskite.params.napfu = 5;
    ca_perovskite.params.molar_mass = 0.1161617;
    ca_perovskite.params.F_0 = -1463358.0;
    ca_perovskite.params.V_0 = 2.745e-05;
    ca_perovskite.params.K_0 = 2.36e11;
    ca_perovskite.params.Kprime_0 = 3.9;
    ca_perovskite.params.G_0 = 1.568315e11;
    ca_perovskite.params.Gprime_0 = 2.22713;
    ca_perovskite.params.debye_0 = 795.779;
    ca_perovskite.params.grueneisen_0 = 1.88839;
    ca_perovskite.params.q_0 = 0.89769;
    ca_perovskite.params.eta_s_0 = 1.28818;
    ca_perovskite.params.equation_of_state = EOSType::SLB3;
  }
};

struct StishoviteFixture {
  Mineral stishovite;
  excesses::ExcessParamVector excess_params_stish;
  StishoviteFixture() {
    stishovite.params.name = "Stishovite";
    stishovite.params.formula = FormulaMap{
          {"Si", 1.0},
          {"O", 2.0}
        };
    stishovite.params.napfu = 3;
    stishovite.params.molar_mass = 0.0600843;
    stishovite.params.F_0 = -818984.6;
    stishovite.params.V_0 = 1.4017e-05;
    stishovite.params.K_0 = 3.143352e11;
    stishovite.params.Kprime_0 = 3.75122;
    stishovite.params.G_0 = 2.2e11;
    stishovite.params.Gprime_0 = 1.93334;
    stishovite.params.debye_0 = 1107.824;
    stishovite.params.grueneisen_0 = 1.37466;
    stishovite.params.q_0 = 2.83517;
    stishovite.params.eta_s_0 = 4.60904;
    stishovite.params.equation_of_state = EOSType::SLB3;
    excess_params_stish = { excesses::LandauParams{-4250.0, 1.0e-9, 0.012} };
    stishovite.set_property_modifier_params(excess_params_stish);
  }
};

struct BridgmaniteSolutionFixture {
  BridgmaniteFixture bdg_fix;
  Solution bdg;
  BridgmaniteSolutionFixture() {
    bdg.set_name("Bridgmanite");
    bdg.set_solution_model(bdg_fix.bdg_solution_model);
    bdg.set_composition(bdg_fix.molar_fractions);
  }
};

struct FerropericlaseSolutionFixture {
  FerropericlaseFixture fp_fix;
  Solution fp;
  FerropericlaseSolutionFixture() {
    fp.set_name("Ferro-periclase");
    fp.set_solution_model(fp_fix.fp_solution_model);
    fp.set_composition(fp_fix.molar_fractions);
  }
};

// Two solutions
struct BdgFperAssemblageFixture {
  BridgmaniteSolutionFixture bdg_fix;
  FerropericlaseSolutionFixture fp_fix;
  Assemblage assemblage;
  BdgFperAssemblageFixture() {
    assemblage.set_name("Bdg + Fper assemblage");
    assemblage.set_averaging_scheme(AveragingType::VRH);
    assemblage.add_phase(bdg_fix.bdg);
    assemblage.add_phase(fp_fix.fp);
    assemblage.set_fractions({0.7, 0.3});
  }
};

// Two minerals
struct CaPvStishAssemblageFixture {
  CaPerovskiteFixture capv_fix;
  StishoviteFixture stish_fix;
  Assemblage assemblage;

  CaPvStishAssemblageFixture() {
    assemblage.set_name("CaPv + Stish assemblage");
    assemblage.set_averaging_scheme(AveragingType::VRH);
    assemblage.add_phase(capv_fix.ca_perovskite);
    assemblage.add_phase(stish_fix.stishovite);
    assemblage.set_fractions({0.6, 0.4});
  }
};

// Two solutions, one mineral
struct PyroliteAssemblageFixture {
  BridgmaniteSolutionFixture bdg_fix;
  FerropericlaseSolutionFixture fp_fix;
  CaPerovskiteFixture capv_fix;
  Assemblage assemblage;

  PyroliteAssemblageFixture() {
    assemblage.set_name("Bdg + Fper + CaPv assemblage");
    assemblage.set_averaging_scheme(AveragingType::VRH);
    assemblage.add_phase(bdg_fix.bdg);
    assemblage.add_phase(fp_fix.fp);
    assemblage.add_phase(capv_fix.ca_perovskite);
    assemblage.set_fractions({0.6, 0.3, 0.1});
  }
};

// One assemblage, one mineral
struct NestedAssemblageFixture {
  PyroliteAssemblageFixture py_fix;
  StishoviteFixture stish_fix;
  Assemblage assemblage;

  NestedAssemblageFixture() {
    assemblage.set_name("Nested assemblage");
    assemblage.set_averaging_scheme(AveragingType::VRH);
    assemblage.add_phase(py_fix.assemblage);
    assemblage.add_phase(stish_fix.stishovite);
    assemblage.set_fractions({0.85, 0.15});
  }
};

#endif // TESTS_SOLUTION_FIXTURES_HPP_INCLUDED
