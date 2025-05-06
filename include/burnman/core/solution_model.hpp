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

#include <string>
#include <map>
#include <Eigen/Dense>
#include "burnman/core/mineral.hpp"

// TODO: rename IdealSolution -> IdealSolutionModel etc.?

/**
 * @class SolutionModel
 * @brief Base class for solution models.
 *
 * Specific solution models should derive from this class and implement
 * all declared functions. Solution models are used by the Solution class.
 */
class SolutionModel {

 public:

  // Using Mineral objects for now - must be instance not derived class!
  // For derived classes we need unique_ptr instead.
  std::vector<Mineral> endmembers;
  /// use endmembers.emplace_back( ) to add

  // Counts
  int n_endmembers;
  int n_sites;
  int n_occupancies;
  // Site multiplicity and occupancy matrices
  Eigen::ArrayXXd site_multiplicities;
  Eigen::ArrayXXd endmember_occupancies;
  Eigen::ArrayXXd endmember_n_occupancies;
  // Chemical formula/site information and strings
  std::vector<std::string> formulas; // Endmember formulas
  std::string empty_formula; // Formula stripped of site info
  std::string general_formula; // Combined solution formula
  std::vector<std::map<std::string, double>> solution_formulae; // Map of site chem for each em.
  std::vector<std::string> site_names; // Generic site names
  std::vector<std::vector<std::string>> sites; // Species on equivalent sites
  // TODO:
  //  Check members and move to protected/private if possible...
  //  Consider size_t for n_endmembers etc.

  virtual ~SolutionModel() = default;
  void process_solution_chemistry();

  // TODO: Size of ArrayXd/MatrixXd in n_endmembers... can we pre-allocate?

  // Public functions always using base class implementation
  /**
   * @brief Compute the excess Gibbs free energy of the solution.
   *
   * @param pressure Pressure at which to evaluate the solution model [Pa].
   * @param temperature Temperature at which to evaluate the solution model [K].
   * @param molar_fractions Molar fractions of the endmembers.
   *
   * @returns Excess Gibbs free energy [J/mol].
   */
  double compute_excess_gibbs_free_energy(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const;

  /**
   * @brief Compute the excess volume of the solution.
   *
   * @param pressure Pressure at which to evaluate the solution model [Pa].
   * @param temperature Temperature at which to evaluate the solution model [K].
   * @param molar_fractions Molar fractions of the endmembers.
   *
   * @returns Excess volume [m^3/mol].
   */
  double compute_excess_volume(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const;

  /**
   * @brief Compute the excess entropy of the solution.
   *
   * @param pressure Pressure at which to evaluate the solution model [Pa].
   * @param temperature Temperature at which to evaluate the solution model [K].
   * @param molar_fractions Molar fractions of the endmembers.
   *
   * @returns Excess entropy [J/K/mol].
   */
  double compute_excess_entropy(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const;

  /**
   * @brief Compute the excess enthalpy of the solution.
   *
   * @param pressure Pressure at which to evaluate the solution model [Pa].
   * @param temperature Temperature at which to evaluate the solution model [K].
   * @param molar_fractions Molar fractions of the endmembers.
   *
   * @returns Excess enthalpy [J/mol].
   */
  double compute_excess_enthalpy(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const;

  // Functions should be made virtual and ovveriden if polynomial
  // solution model added. All three return 0 at present.
  /**
   * @brief Compute excess heat capacity at current state.
   *
   * @returns Excess heat capacity [J/K/mol].
   */
  double compute_Cp_excess() const;

  /**
   * @brief Compute excess alpha*V at current state.
   *
   * @returns Excess in [m^3/K/mol].
   */
  double compute_alphaV_excess() const;

  /**
   * @brief Compute excess V/K_T at current state.
   *
   * @returns Excess in [m^3/Pa/mol].
   */
  double compute_VoverKT_excess() const;

  // Pure virtual functions to be ovveriden in derived classes
  /**
   * @brief Compute the excess Gibbs free energy for each endmember.
   *
   * @param pressure Pressure at which to evaluate the solution model [Pa].
   * @param temperature Temperature at which to evaluate the solution model [K].
   * @param molar_fractions Molar fractions of the endmembers.
   *
   * @returns Excess partial Gibbs free energies [J/mol].
   */
  virtual Eigen::ArrayXd compute_excess_partial_gibbs_free_energies(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const = 0;

  /**
   * @brief Compute the excess entropy for each endmember.
   *
   * @param pressure Pressure at which to evaluate the solution model [Pa].
   * @param temperature Temperature at which to evaluate the solution model [K].
   * @param molar_fractions Molar fractions of the endmembers.
   *
   * @returns Excess partial entropies [J/K/mol].
   */
  virtual Eigen::ArrayXd compute_excess_partial_entropies(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const = 0;

  /**
   * @brief Compute the excess partial volume for each endmember.
   *
   * @param pressure Pressure at which to evaluate the solution model [Pa].
   * @param temperature Temperature at which to evaluate the solution model [K].
   * @param molar_fractions Molar fractions of the endmembers.
   *
   * @returns Excess partial volumes [m^3/mol].
   */
  virtual Eigen::ArrayXd compute_excess_partial_volumes(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const = 0;

  /**
   * @brief Compute the activities of each endmember.
   *
   * @param pressure Pressure at which to evaluate the solution model [Pa].
   * @param temperature Temperature at which to evaluate the solution model [K].
   * @param molar_fractions Molar fractions of the endmembers.
   *
   * @returns Activities [dimensionless].
   */
  virtual Eigen::ArrayXd compute_activities(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const = 0;

  /**
   * @brief Compute the activity coefficients of the endmembers.
   *
   * @param pressure Pressure at which to evaluate the solution model [Pa].
   * @param temperature Temperature at which to evaluate the solution model [K].
   * @param molar_fractions Molar fractions of the endmembers.
   *
   * @returns Activity coefficients [dimensionless].
   */
  virtual Eigen::ArrayXd compute_activity_coefficients(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const = 0;

  /**
   * @brief Compute second compositional derivative of the Gibbs free energy.
   *
   * @param pressure Pressure at which to evaluate the solution model [Pa].
   * @param temperature Temperature at which to evaluate the solution model [K].
   * @param molar_fractions Molar fractions of the endmembers.
   *
   * @returns Hessian of the Gibbs free energy [J].
   */
  virtual Eigen::MatrixXd compute_gibbs_hessian(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const = 0;

  /**
   * @brief Compute second compositional derivative of the entropy.
   *
   * @param pressure Pressure at which to evaluate the solution model [Pa].
   * @param temperature Temperature at which to evaluate the solution model [K].
   * @param molar_fractions Molar fractions of the endmembers.
   *
   * @returns Hessian of the entropy [J/K].
   */
  virtual Eigen::MatrixXd compute_entropy_hessian(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const = 0;

  /**
   * @brief Compute second compositional derivative of the volume.
   *
   * @param pressure Pressure at which to evaluate the solution model [Pa].
   * @param temperature Temperature at which to evaluate the solution model [K].
   * @param molar_fractions Molar fractions of the endmembers.
   *
   * @return Hessian of the partial volumes [m^3].
   */
  virtual Eigen::MatrixXd compute_volume_hessian(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const = 0;

 protected:
  ;
 private:
  ;
};

/**
 * @class IdealSolution
 * @brief Ideal solution model.
 *
 * Derived from SolutionModel.
 * Calculates the excess gibbs free energy and etropy due to configurational entropy.
 * Excess internal energy and volume are equal to zero.
 *
 * The multiplicity of each type of site in the structure is allowed to change
 * linearly as a function of endmember proportions. This makes the model equivalent
 * to the entropic part of a Temkin-type model (Temkin, 1945).
 *
 */
class IdealSolution : public SolutionModel{

 public:

  // Init - process_solution_chemistry
  //      - calc endmember_configurational_entropy

  // Public functions overriden from base class
  Eigen::ArrayXd compute_excess_partial_gibbs_free_energies(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::ArrayXd compute_excess_partial_entropies(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::ArrayXd compute_excess_partial_volumes(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::ArrayXd compute_activities(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::ArrayXd compute_activity_coefficients(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::MatrixXd compute_gibbs_hessian(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::MatrixXd compute_entropy_hessian(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::MatrixXd compute_volume_hessian(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

 protected:

 private:

  // Public or protected member variable?
  Eigen::ArrayXd endmember_configurational_entropies;
  Eigen::ArrayXd compute_endmember_configurational_entropies() const;

  // This is unused in python?
  double compute_configurational_entropy(
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::ArrayXd compute_ideal_excess_partial_gibbs(
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::ArrayXd compute_ideal_excess_partial_entropies(
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::ArrayXd compute_ideal_activities(
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::ArrayXd compute_log_ideal_activities(
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::MatrixXd compute_log_ideal_activity_derivatives(
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::MatrixXd compute_ideal_entropy_hessian(
    const Eigen::ArrayXd& molar_fractions) const;

  // Want to make and cache, ones, eyeones, eye
  // Not really used --> would make protected and create in setup

};

/**
 * @class AsymmetricRegularSolution
 * @brief Asymmetric regular solution model.
 *
 * Derived from IdealSolution.
 * Implements the asymmetric regular solution model desribed in Holland &
 * Powell, 2003.
 *
 * The excess non-configurational Gibbs energy is given by the expression:
 * \f[
 * \mathcal{G}_{\textrm{excess}} = \alpha^T p (\phi^T W \phi)
 * \f]
 * \f$\alpha\f$ is a vector of van Laar parameters governing asymmetry in the
 * excess properties, i.e.:
 * \f[
 * \phi_i = \frac{\alpha_i p_i}{\sum_{k=1}^{n} \alpha_k p_k}, \quad
 * W_{ij} = \frac{2 w_{ij}}{\alpha_i + \alpha_j} \quad \text{for } i < j
 * \f]
 *
 * See `SymmetricRegularSolution' for the special case where all alpha = 1.
 */
class AsymmetricRegularSolution : public IdealSolution{

 public:

  // Setup in __init__!

  // Alphas (van Laar) and interaction parameters
  Eigen::ArrayXd alphas;
  Eigen::MatrixXd W_e; // energy interactions
  Eigen::MatrixXd W_s; // entropy interactions
  Eigen::MatrixXd W_v; // volume interactions

  Eigen::ArrayXd compute_excess_partial_gibbs_free_energies(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::ArrayXd compute_excess_partial_entropies(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::ArrayXd compute_excess_partial_volumes(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::ArrayXd compute_activities(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::ArrayXd compute_activity_coefficients(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::MatrixXd compute_gibbs_hessian(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::MatrixXd compute_entropy_hessian(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

  Eigen::MatrixXd compute_volume_hessian(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const override;

 protected:

 private:

  Eigen::ArrayXd compute_phi(
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::ArrayXd compute_non_ideal_interactions(
    const Eigen::MatrixXd& W,
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::ArrayXd compute_non_ideal_excess_partial_gibbs(
    double pressure,
    double temperature,
    const Eigen::ArrayXd& molar_fractions) const;

  Eigen::MatrixXd compute_non_ideal_hessian(
    const Eigen::MatrixXd& interactions,
    const Eigen::ArrayXd& molar_fractions) const;

};

/**
 * @class SymmetricRegularSolution
 * @brief Convenience class for a symmetric regular solution model.
 *
 * This class is a special case of `AsymmetricRegularSolution' with all
 * alphas set to 1.
 *
 */
class SymmetricRegularSolution : public AsymmetricRegularSolution{
  // Only change to setup, with alpha = 1,1,1,...
  ;
};

#endif // BURNMAN_CORE_SOLUTION_MODEL_HPP_INCLUDED

// empty_formula [string] --> abbreviated formula with empty [] for sites
// general_formula --> formula with comma separated possible species
// n_sites [int]
// sites --> list of list of strings (check best way to implement, may ignore if unused)
// site_names --> vector of strings
// n_occupancies [int] --> sum of number of possible species
//  e.g. [[A][B], [B][C1/2D1/2]] would = 5 (2 on 1, 3 on 2)
// site_multiplicities [2D array] --> 1D for each endmember
// endmember_occupancies [2D array] --> as above
// endmember_noccupancies [2D array] --> number of atoms instead of fraction

// TODO:
//  Base class, ideal, asymmetric
//  Public getters and protected compute functions
//  Make virtual where needed to override in derived classes
/*
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