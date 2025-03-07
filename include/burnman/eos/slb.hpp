/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_EOS_SLB_HPP_INCLUDED
#define BURNMAN_EOS_SLB_HPP_INCLUDED

#include "burnman/core/equation_of_state.hpp"
#include "burnman/utils/eos.hpp"

/**
 * @class SLB3
 * @brief Class for the Mie-Grueneisen-Debye equation of state detailed
 *        in Stixrude and Lithgow-Bertelloni, 2005 (Geophys. J. Int.).
 *
 * Uses third order finite strain expansion for the shear modulus.
 * To fit shear modulus to second order expansion @see SLB2.
 *
 * @note All functions assume SI units for all properties.
 */
class SLB3 : public EquationOfState{
 public:

  // Helper functions
  bool validate_parameters(MineralParams& params) override;

  // Specific EOS functions
  double compute_volume(
    double pressure,
    double temperature,
    const MineralParams& params) const override;

  double compute_pressure(
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_grueneisen_parameter(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_isothermal_bulk_modulus_reuss(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_isentropic_bulk_modulus_reuss(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_shear_modulus(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_molar_heat_capacity_v(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_molar_heat_capacity_p(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_thermal_expansivity(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_gibbs_free_energy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_entropy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_molar_internal_energy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_helmholtz_free_energy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

  double compute_enthalpy(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;

 protected:

  /**
   * @brief Helper function for shear modulus calculation
   *
   * @param temperatue Temp in [K].
   * @param volume Vol in [m^3].
   * @param params Mineral parameters object.
   *
   * @return delta_G portion of G
   */
   double compute_shear_modulus_delta(
    double temperatue,
    double volume,
    const MineralParams& params) const;

 private:

  /**
   * @brief Computes the volume dependent Grueneisen parameter.
   *
   * @param x V_0/V.
   * @param params Mineral parameters object of type MineralParams
   * 
   * @return Grueneisen parameter [unitless].
   */
  static double compute_slb_grueneisen_parameter(
    double x,
    const MineralParams& params);

  /**
   * @brief Compute the Debye temperature.
   *
   * Finite strain approximation.
   *
   * @param x V_0/V (inverse compression).
   * @param params Mineral parameters object of type MineralParams.
   *
   * @return Debye temperature in [K].
   */
  static double compute_debye_temperature(
    double x,
    const MineralParams& params);

  /**
   * @brief Compute q
   *
   * Finite strain approximation for the isotropic volume strain
   * derivative of the grueneisen parameter.
   *
   * @param x V_0/V (inverse compression).
   * @param params Mineral parameters object.
   *
   * @return q [unitless]
   */
  static double compute_volume_dependent_q(
    double x,
    const MineralParams& params);

  /**
   * @brief Compute eta_s0
   *
   * Finite strain approximation for the isotropic shear strain
   * derivative of the grueneisen parameter.
   *
   * @param x V_0/V (inverse compression).
   * @param params Mineral parameters object.
   *
   * @return eta_s [unitless]
   */
  static double compute_isotropic_eta_s(
    double x,
    const MineralParams& params);

  /**
   * @brief GSL function wrapper to compute P(V) - P
   * 
   * @param x Volume to test (passed by solver)
   * @param p Generic pointer for parameter object
   * @see `ParamsGSL::SolverParams_SLB`
   */
  static double slb_gsl_wrapper(double x, void* p);

  // TODO: wrapper for K_T for GSL Brent also

};


/**
 * @class SLB2
 * @brief SLB EOS with second order expansion of shear modulus.
 *
 * Derived from SLB3.
 * Overrides compute_shear_modulus to use the second order expansion.
 * Use with caution - prefer SLB3.
 *
 * @note All functions assume SI units for all properties.
 */
class SLB2 : public SLB3{
 public:
  /**
   * @copydoc EquationOfState::compute_shear_modulus
   *
   * @note Second order expansion
   */
  double compute_shear_modulus(
    double pressure,
    double temperature,
    double volume,
    const MineralParams& params) const override;
};

/**
 * @class SLB3Conductive
 * @brief SLB EOS with third order expansion of shear modulus and a
 *        contribution to the Helmholtz free energy that arises from the
 *        thermal excitation of electrons (Bukowinski, 1977).
 *
 * Derived from SLB3.
 * Overrides -->
 *
 * @note All functions assume SI units for all properties.
 */
class SLB3Conductive : public SLB3{
 public:

  // Functions to override:
  //  compute_volume
  //  compute_pressure
  //  compute_isothermal_bulk_modulus_reuss
  //  compute_molar_heat_capacity_v
  //  compute_thermal_expansivity
  //  compute_entropy
  //  compute_helmholtz_free_energy
  //  compute_grueneisen_parameter
  //  validate_parameters

};


#endif // #define BURNMAN_EOS_SLB_HPP_INCLUDED



