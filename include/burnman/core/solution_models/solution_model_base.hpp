/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_SOLUTION_MODELS_BASE_HPP_INCLUDED
#define BURNMAN_CORE_SOLUTION_MODELS_BASE_HPP_INCLUDED

#include <map>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/core/mineral.hpp"

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
  // Using Eigen::Index (usually std::ptrdiff_t) - cast to size_t for STL containers
  Eigen::Index n_endmembers;
  Eigen::Index n_sites;
  Eigen::Index n_occupancies;
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
  // const std::vector<Mineral>& minerals() const { return minerals; } ...etc.
  //  Consider size_t for n_endmembers etc.

  // Constructor
  SolutionModel(const PairedEndmemberList& endmember_list);

  virtual ~SolutionModel() = default;
  void process_solution_chemistry();

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

#endif // BURNMAN_CORE_SOLUTION_MODEL_BASE_HPP_INCLUDED
