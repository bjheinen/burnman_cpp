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
//#include "burnman/core/solution.hpp"
#include "burnman/core/mineral.hpp"
#include "burnman/core/solution_model.hpp"
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

#endif // TESTS_SOLUTION_FIXTURES_HPP_INCLUDED
