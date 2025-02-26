/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_SOLUTION_MODEL_HPP_INCLUDED
#define BURNMAN_CORE_SOLUTION_MODEL_HPP_INCLUDED

/**
 * Base class for solution models.

  TODO
  
 */ 
class SolutionModel {

 public:

  virtual ~SolutionModel() = default;

 protected:
  ;
 private:
  ;

};

#endif // BURNMAN_CORE_SOLUTION_MODEL_HPP_INCLUDED

/*
// Derived classes
--> TODO: AsymmetricRegularSolution
MechanicalSolution
IdealSolution
AsymmetricRegularSolution
SymmetricRegularSolution --> special case of asym with alphas = ones
SubRegularSolution
FunctionSolution
PolynomialSolution

// Public getters and protected compute functions 
// Make virtual where needed to override in derived classes
excess_gibbs_free_energy
excess_partial_gibbs_free_energies
excess_volume
excess_partial_volumes
excess_entropy
excess_partial_entropies
excess_enthalpy
Cp_excess
alphaV_excess
VoverKT_excess
// In derived
activity_coefficients
activities
// Asymmetric
Constructor logic
phi
non_ideal_interactions
non_ideal_excess_partial_gibbs
excess_partial_gibbs_free_energies --> override
excess_partial_entropies --> override
excess_partial_volumes --> override
gibbs_hessian
entorpy_hessian
volume_hessian
*/
